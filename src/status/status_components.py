import dash_bootstrap_components as dbc
from dash import dcc, html

def status_progress():
	m = dbc.Progress(id="status_progress")
	return m

def status_history():
	m = dbc.Row(children=[
			dbc.Col(children=[
				html.Div(id="status_history", children=[""])
			]),
		])
	return m

def status_state():
	m = dbc.Row(children=[
			dbc.Col(children=[
				html.Div(id="status_state", children=[""])
			]),
		])
	return m
