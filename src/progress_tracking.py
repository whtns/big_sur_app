"""
Progress tracking module using smtplib for sending email notifications.
Integrates with Celery tasks to track and report progress.
"""

import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from datetime import datetime
import json
import redis
import os
from typing import Optional, Dict, Any

from layouts import demo
if (demo == False):
    try:
        from configmodule.production_config import ProductionConfig as FlaskConfig
    except:
        from configmodule.default_config import DefaultConfig as FlaskConfig
else:
    from configmodule.default_config import DefaultConfig as FlaskConfig


class ProgressTracker:
    """
    Tracks task progress and sends email updates via smtplib.
    Stores progress state in Redis for real-time access.
    """
    
    def __init__(self):
        self.config = FlaskConfig()
        # respect environment / compose settings; default to redis service name
        redis_host = os.environ.get('REDIS_HOST', os.environ.get('CELERY_BROKER_HOST', 'redis'))
        redis_port = int(os.environ.get('REDIS_PORT', os.environ.get('CELERY_BROKER_PORT', 6379)))
        redis_db = int(os.environ.get('REDIS_DB', 0))
        self.redis_client = redis.Redis(host=redis_host, port=redis_port, db=redis_db)
        self.smtp_server = getattr(self.config, 'SMTP_SERVER', 'localhost')
        self.smtp_port = getattr(self.config, 'SMTP_PORT', 587)
        self.smtp_username = getattr(self.config, 'SMTP_USERNAME', None)
        self.smtp_password = getattr(self.config, 'SMTP_PASSWORD', None)
        self.smtp_sender = getattr(self.config, 'SMTP_SENDER', 'noreply@bigsur.local')
        self.use_tls = getattr(self.config, 'SMTP_USE_TLS', True)
    
    def update_progress(self, task_id: str, user_email: str, 
                       task_name: str, progress: int, 
                       total: int, details: Optional[str] = None) -> None:
        """
        Update task progress and store in Redis.
        
        Args:
            task_id: Unique task identifier
            user_email: Email of user running the task
            task_name: Name of the task being tracked
            progress: Current progress count
            total: Total items to process
            details: Optional details about current step
        """
        progress_key = f"progress:{task_id}"
        progress_data = {
            'task_id': task_id,
            'task_name': task_name,
            'user_email': user_email,
            'progress': progress,
            'total': total,
            'percentage': int((progress / total * 100)) if total > 0 else 0,
            'timestamp': datetime.now().isoformat(),
            'details': details or ''
        }
        
        # Store in Redis with 1 hour expiry
        self.redis_client.setex(
            progress_key, 
            3600, 
            json.dumps(progress_data)
        )
    
    def get_progress(self, task_id: str) -> Optional[Dict[str, Any]]:
        """
        Retrieve progress for a specific task.
        
        Args:
            task_id: Unique task identifier
            
        Returns:
            Progress data dictionary or None if not found
        """
        progress_key = f"progress:{task_id}"
        data = self.redis_client.get(progress_key)
        if data:
            return json.loads(data)
        return None
    
    def send_progress_email(self, recipient_email: str, task_name: str,
                           progress: int, total: int, 
                           details: Optional[str] = None,
                           task_status: str = 'in_progress') -> bool:
        """
        Send progress update email via SMTP.
        
        Args:
            recipient_email: Email address to send to
            task_name: Name of the task
            progress: Current progress count
            total: Total items to process
            details: Optional details
            task_status: Status of task ('in_progress', 'completed', 'error')
            
        Returns:
            True if email sent successfully, False otherwise
        """
        try:
            percentage = int((progress / total * 100)) if total > 0 else 0
            
            # Create message
            msg = MIMEMultipart('alternative')
            msg['Subject'] = f"[BigSuR] {task_name} Progress Update - {percentage}%"
            msg['From'] = self.smtp_sender
            msg['To'] = recipient_email
            
            # Status emoji/text
            status_emoji = {
                'in_progress': '⏳',
                'completed': '✅',
                'error': '❌'
            }
            
            # Plain text version
            text = f"""
BigSuR Task Progress Update

Task: {task_name}
Status: {task_status.replace('_', ' ').title()}
Progress: {progress}/{total} ({percentage}%)
Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

{f'Details: {details}' if details else ''}

This is an automated notification from BigSuR.
            """.strip()
            
            # HTML version
            html = f"""
<html>
  <body>
    <div style="font-family: Arial, sans-serif; max-width: 600px; margin: 0 auto;">
      <h2 style="color: #333;">BigSuR Task Progress Update</h2>
      <div style="background-color: #f5f5f5; padding: 15px; border-radius: 5px; margin: 10px 0;">
        <p><strong>Task:</strong> {task_name}</p>
        <p><strong>Status:</strong> {status_emoji.get(task_status, '⏳')} {task_status.replace('_', ' ').title()}</p>
        <div style="margin: 15px 0;">
          <strong>Progress: {progress}/{total} ({percentage}%)</strong>
          <div style="background-color: #ddd; border-radius: 5px; height: 20px; overflow: hidden;">
            <div style="background-color: #4CAF50; height: 100%; width: {percentage}%; transition: width 0.3s;"></div>
          </div>
        </div>
        <p><strong>Timestamp:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        {f'<p><strong>Details:</strong> {details}</p>' if details else ''}
      </div>
      <p style="color: #888; font-size: 12px;">This is an automated notification from BigSuR.</p>
    </div>
  </body>
</html>
            """.strip()
            
            part1 = MIMEText(text, 'plain')
            part2 = MIMEText(html, 'html')
            msg.attach(part1)
            msg.attach(part2)
            
            # Send email
            if self.use_tls:
                server = smtplib.SMTP(self.smtp_server, self.smtp_port)
                server.starttls()
                if self.smtp_username and self.smtp_password:
                    server.login(self.smtp_username, self.smtp_password)
            else:
                server = smtplib.SMTP_SSL(self.smtp_server, self.smtp_port)
                if self.smtp_username and self.smtp_password:
                    server.login(self.smtp_username, self.smtp_password)
            
            server.sendmail(self.smtp_sender, recipient_email, msg.as_string())
            server.quit()
            
            print(f"[INFO] Progress email sent to {recipient_email}")
            return True
            
        except Exception as e:
            print(f"[ERROR] Failed to send progress email: {str(e)}")
            return False
    
    def send_completion_email(self, recipient_email: str, task_name: str,
                             success: bool = True, 
                             result_details: Optional[str] = None) -> bool:
        """
        Send task completion email via SMTP.
        
        Args:
            recipient_email: Email address to send to
            task_name: Name of the completed task
            success: Whether task completed successfully
            result_details: Optional details about the result
            
        Returns:
            True if email sent successfully, False otherwise
        """
        try:
            status = "Completed" if success else "Failed"
            status_emoji = "✅" if success else "❌"
            
            msg = MIMEMultipart('alternative')
            msg['Subject'] = f"[BigSuR] {task_name} {status} {status_emoji}"
            msg['From'] = self.smtp_sender
            msg['To'] = recipient_email
            
            # Plain text version
            text = f"""
BigSuR Task Completion Notification

Task: {task_name}
Status: {status}
Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

{f'Details: {result_details}' if result_details else ''}

This is an automated notification from BigSuR.
            """.strip()
            
            # HTML version
            html = f"""
<html>
  <body>
    <div style="font-family: Arial, sans-serif; max-width: 600px; margin: 0 auto;">
      <h2 style="color: #333;">BigSuR Task Completion</h2>
      <div style="background-color: {'#d4edda' if success else '#f8d7da'}; padding: 15px; border-radius: 5px; margin: 10px 0; border-left: 4px solid {'#28a745' if success else '#dc3545'};">
        <p><strong>Task:</strong> {task_name}</p>
        <p><strong>Status:</strong> {status_emoji} <span style="color: {'#28a745' if success else '#dc3545'}; font-weight: bold;">{status}</span></p>
        <p><strong>Timestamp:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        {f'<p><strong>Details:</strong> {result_details}</p>' if result_details else ''}
      </div>
      <p style="color: #888; font-size: 12px;">This is an automated notification from BigSuR.</p>
    </div>
  </body>
</html>
            """.strip()
            
            part1 = MIMEText(text, 'plain')
            part2 = MIMEText(html, 'html')
            msg.attach(part1)
            msg.attach(part2)
            
            # Send email
            if self.use_tls:
                server = smtplib.SMTP(self.smtp_server, self.smtp_port)
                server.starttls()
                if self.smtp_username and self.smtp_password:
                    server.login(self.smtp_username, self.smtp_password)
            else:
                server = smtplib.SMTP_SSL(self.smtp_server, self.smtp_port)
                if self.smtp_username and self.smtp_password:
                    server.login(self.smtp_username, self.smtp_password)
            
            server.sendmail(self.smtp_sender, recipient_email, msg.as_string())
            server.quit()
            
            print(f"[INFO] Completion email sent to {recipient_email}")
            return True
            
        except Exception as e:
            print(f"[ERROR] Failed to send completion email: {str(e)}")
            return False


# Create a singleton instance
progress_tracker = ProgressTracker()
