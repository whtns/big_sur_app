import os
class DefaultConfig(object):
	SECRET_KEY = os.environ.get("SECRET_KEY", b"\xc9j\xa2@k\x04\x0e\x8a\xe9\xb6\xfbA\xdfsU\x05\xdfe\xec@\x05\x0b\xfd\x9a")

	SECURITY_PASSWORD_SALT = '1381971702542417111521178241349271342290814913196'
	SECURITY_REGISTERABLE = True #allows user registration
	DATABASE_URI = 'sqlite:////app/data/misc.db' # user security database
	
	# SMTP Configuration for progress tracking emails
	SMTP_SERVER = os.environ.get('SMTP_SERVER', 'localhost')
	SMTP_PORT = int(os.environ.get('SMTP_PORT', 587))
	SMTP_USERNAME = os.environ.get('SMTP_USERNAME', None)
	SMTP_PASSWORD = os.environ.get('SMTP_PASSWORD', None)
	SMTP_SENDER = os.environ.get('SMTP_SENDER', 'bigsur@localhost')
	SMTP_USE_TLS = os.environ.get('SMTP_USE_TLS', 'True').lower() == 'true'