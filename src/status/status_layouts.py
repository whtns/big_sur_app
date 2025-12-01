import dash
import dash_bootstrap_components as dbc
from dash import dcc, html

from . import status_components as cc

def status_layout():
	m = dbc.Card(children=[
			dbc.CardHeader(cc.status_progress()),
			dbc.CardBody(cc.status_history()),
			dbc.CardFooter(cc.status_state())
		])
	return m