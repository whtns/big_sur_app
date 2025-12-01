import dash
from dash import dcc, html
import dash_bootstrap_components as dbc

from . import importing_components as cc

def importing_layout(demo=False):
	m = dbc.Container(children=[
			dbc.Row(children=[
				dbc.Col(children=[
					cc.importing_greeting(demo=demo)
				])
			]),
			dbc.Row(children=[
		    	dbc.Col(children=[
					# CardDeck is deprecated in newer dbc versions — use grid of Columns
					dbc.Row(children=[
						dbc.Col(cc.importing_dataset_dropdown(), width=4),
						dbc.Col(cc.importing_user_dataset_list(demo=demo), width=4),
						dbc.Col(cc.importing_data_upload(demo=demo), width=4),
					]),
			    ])
			], id="upload-collapse")
		], fluid=True)
	return m