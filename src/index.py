import os
# Ensure Numba JIT is disabled early to avoid llvmlite JIT/runtime
# initialization before Gunicorn forks workers.
os.environ.setdefault('NUMBA_DISABLE_JIT', '1')

import dash
from dash import dcc, html
from dash.dependencies import Input, Output

import flask
import flask_security
from flask_security import current_user, logout_user

import uuid

from app import app
import layouts

app.layout = layouts.main_layout

from inputoutput import inputoutput_callbacks
from processing import processing_callbacks
from markergenes import markergenes_callbacks
from annotation import annotation_callbacks
from importing import importing_callbacks
from exporting import exporting_callbacks
from status import status_callbacks
from hvg import hvg_callbacks
from correlations.correlation_callbacks import register_correlation_callbacks

# Register correlation callbacks
register_correlation_callbacks(app)

server = app.server

# extra static routes
@server.route('/favicon.ico')
def send_favicon():
    return flask.send_from_directory("/Library/WebServer/Documents/BigSuR/assets/", 
    								 "favicon.ico")

if __name__ == '__main__':
    # Dash 2.10+ replaced `run_server` with `run` (ObsoleteAttributeException otherwise)
    # Keep the same debug and host arguments so running `python src/index.py` works.
    try:
        app.run(debug=True, host='0.0.0.0')
    except AttributeError:
        # Fallback for older Dash versions that still use run_server
        app.run_server(debug=True, host='0.0.0.0')