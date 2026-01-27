# celery configuration: prefer environment variables set by compose, with sensible defaults
import os
from celery import Celery

broker_url = os.environ.get('CELERY_BROKER_URL', 'redis://redis:6379/0')
backend_url = os.environ.get('CELERY_RESULT_BACKEND', broker_url)

task_queue = Celery('BigSuR', broker=broker_url, backend=backend_url, include=['tasks.tasks'])
# Optional configuration, see the application user guide.
task_queue.conf.update(
    result_expires=3600,
)
if __name__ == '__main__':
    task_queue.start()