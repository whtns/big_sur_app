import dash
import ipdb

from flask import Flask
from flask_mail import Mail
from flask_caching import Cache
from flask_security import Security, login_required, \
     SQLAlchemySessionUserDatastore
from flask_admin import Admin
from sqlalchemy import exc
from flask_bootstrap import Bootstrap4

from user_management.database import db_session, init_db
from user_management.models import User, Role

from admin.views import AdminView

import layouts
from tasks.tasks import send_flask_mail
from layouts import demo
if (demo is False):
    try:
    	from configmodule.production_config import ProductionConfig as FlaskConfig
    except:
    	from configmodule.default_config import DefaultConfig as FlaskConfig
else:
	from configmodule.default_config import DefaultConfig as FlaskConfig

security = Security()

# instantiate the app with static file paths
import os
static_path = os.path.join(os.path.dirname(__file__), 'assets')
server = Flask(__name__, static_url_path='/assets', static_folder=static_path)
server.config.from_object(FlaskConfig())
Bootstrap4(server)

app = dash.Dash(server=server, show_undo_redo=False,
                url_base_pathname="/", update_title=None)
app.title = "BigSuR"
app.config.suppress_callback_exceptions = True

# Setup cache
import os

# Configure Redis URL for Flask-Caching from environment, default to compose service
redis_url = os.environ.get('CACHE_REDIS_URL', os.environ.get('REDIS_URL', 'redis://redis:6379/0'))
cache = Cache(app.server, config={
	'CACHE_TYPE': 'redis',
	'CACHE_REDIS_URL': redis_url,
	'CACHE_THRESHOLD': 200,  # should be set to max number of active users
	"CACHE_DEFAULT_TIMEOUT": 30000
	}
)

# Setup Flask-Security
user_datastore = SQLAlchemySessionUserDatastore(db_session,
                                                User, Role)
security_ctx = security.init_app(app.server, user_datastore, register_blueprint=True)


def delay_flask_security_mail(msg):
	print("[DEBUG] running delay_flask_security_mail")
	with app.server.app_context():
		send_flask_mail.delay(subject=msg.subject, sender=msg.sender,
		    	              recipients=msg.recipients, body=msg.body,
		    	              html=msg.html)
	return None
try:
	# security.init_app may return None depending on Flask-Security version; use the extension instance
	# to register the send_mail_task callback to avoid AttributeError
	target = security_ctx if security_ctx is not None else security
	target.send_mail_task(delay_flask_security_mail)
except AttributeError:
	# If the extension API differs, silently continue (mail will be sent synchronously as fallback)
	pass

# Create a user to test with
def initialize_database():
	init_db()
	db_session.commit()

# Register the initialization function in a way that's compatible with
# Flask 1/2 (`before_first_request`) and Flask 3+ (`before_serving`). If
# neither hook exists, run initialization immediately as a fallback.
if hasattr(server, 'before_first_request'):
	server.before_first_request(initialize_database)
elif hasattr(server, 'before_serving'):
	server.before_serving(initialize_database)
else:
	# Last-resort: call immediately (useful for some test runners)
	initialize_database()
'''
def create_test_user():
    init_db()
    try:
    	user_datastore.create_user(email='nigeil@yahoo.com', username="nigeil", password='01010101')
    except exc.SQLAlchemyError:
    	print("[DEBUG] user already exists")
    db_session.commit()
'''

# Setup the admin interface
admin = Admin(app.server)
admin.add_view(AdminView(User, db_session, endpoint="user_admin"))
admin.add_view(AdminView(Role, db_session, endpoint="role_admin"))