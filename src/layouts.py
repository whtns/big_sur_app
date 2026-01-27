from dash import html, dcc
import dash_bootstrap_components as dbc

from flask_security import current_user

import uuid
import base64
import os

from processing.processing_layouts import *
from markergenes.markergenes_layouts import *
from annotation.annotation_layouts import *
from importing.importing_layouts import *
from exporting.exporting_layouts import *
from status.status_layouts import *
from hvg.hvg_layouts import *
from correlations.correlation_layouts import correlation_tab_layout

demo = False 

def titlebar_layout(demo=False):	
	# load the logo
	# Try multiple paths: production path, development path, assets fallback
	logo_paths = [
		"/Library/WebServer/Documents/BigSuR/assets/img/BigSuR_logo.png",
		"./images/BigSuR_logo.png",
		"/app/images/BigSuR_logo.png"  # Docker path
	]
	
	logo_src = None
	for logo_path in logo_paths:
		if os.path.exists(logo_path):
			try:
				with open(logo_path, 'rb') as _f:
					BigSuR_logo = base64.b64encode(_f.read()).decode()
				logo_src = f"data:image/png;base64,{BigSuR_logo}"
				break
			except Exception:
				continue
	
	# Fall back to assets/images path if file not found
	if logo_src is None:
		logo_src = "/assets/BigSuR_logo.png"

	if (demo == True):
		login_logout = html.Span()
	elif (current_user and current_user.is_authenticated):
		login_logout = html.A("Logout", href="/logout", className="btn btn-dark text-light", role="button")
	else:
		login_logout = html.A("Login", href="/login", className="btn btn-dark text-light", role="button")

	# Build navbar using Dash/Bootstrap components instead of raw HTML string
	return html.Nav(
		className="navbar navbar-light navbar-expand-md bg-light",
		style={"fontFamily": "'Open Sans', sans-serif"},
		children=[
			html.Div(
				className="container-fluid",
				children=[
						html.A(html.Img(src=logo_src, className="img-fluid", width="200", height="200"), href="/about"),
					html.A(className="navbar-brand", href="#"),
					html.Button(html.Span(className="navbar-toggler-icon"), **{"data-toggle": "collapse", "data-target": "#navcol-1", "className": "navbar-toggler"}),
					html.Div(
						className="collapse navbar-collapse",
						id="navcol-1",
						children=[
							html.Ul([
								html.Li(html.A("About", href="/about", className="nav-link active"), className="nav-item"),
								html.Li("", className="nav-item"),
								html.Li(html.A("Docs", href="/docs", className="nav-link"), className="nav-item"),
								html.Li(html.A("Code", href="/code", className="nav-link"), className="nav-item")
							], className="nav navbar-nav"),
							html.Ul([html.Li(login_logout, className="nav-item")], className="nav navbar-nav ml-auto")
						]
					)
				]
			)
		]
	)

'''
	ret = html.Nav(children=[
			html.A(children=[
				html.Img(src="data:image/png;base64,{}".format(BigSuR_logo.decode()), height="45px"),
				html.H3("A Multi-informatic Cellular Visualization tool"),
			], href="https://BigSuR.works", className="navbar-brand"),
			html.Div(children=[
				dbc.Nav(children=[
					(html.A("Logout", href="/logout", className="nav-item nav-link") if (demo is False) else html.A("Login", href="/login", className="nav-item nav-link")),
					html.A("About", href="/about", className="nav-item nav-link"),
					html.A("Documentation", href="/documentation", className="nav-item nav-link")
			    ]),
			], className="navbar-expand flex-grow-1 text-left", id="navbarNav")
	], className="navbar navbar-light bg-light")

'''

def main_layout():
	if ((current_user) and (current_user.is_authenticated is True)):
		session_id = str(current_user.id)
	else:
		session_id = str(uuid.uuid4())

	return html.Div(children=[
	    html.Div(session_id, id='session-id', style={'display': 'none'}),
	    html.Div([], id='null_container_0', style={'display': 'none'}), # stores nothing, useful output endpoint for callbacks that don't modify anything in the UI
	    html.Div([], id='null_container_1', style={'display': 'none'}), # stores nothing, useful output endpoint for callbacks that don't modify anything in the UI
	    html.Div([], id='null_container_2', style={'display': 'none'}), # stores nothing, useful output endpoint for callbacks that don't modify anything in the UI

	    dcc.Interval(id="cleanup_interval",
	    			 interval=(15 * 60 * 1000),
	    			 max_intervals=-1),

		dcc.Interval(id="status_interval", 
					 interval=(1 * 2 * 1000),
					 max_intervals=-1),

		# title matter
		titlebar_layout(demo=demo),

	    dbc.Container(fluid=True, children=[
		    dbc.Row(children=[
		    	# Main layout
		    	dbc.Col(children=[
		    		dbc.Tabs(children=[
		    			### IMPORTING TAB ###
			        	dbc.Tab(importing_layout(demo=demo), label="Load", 
			        			tab_id="importing_tab"),

				    	### PROCESSING TAB ###
			        	dbc.Tab(processing_layout(demo=demo), label="Preprocess", tab_id="processing_tab"),
						
						### MARKER GENE TAB ###
						dbc.Tab(markergenes_layout(), label="Marker genes", tab_id="markergenes_tab"),
					
						### HIGHLY VARIABLE GENES TAB ###
						dbc.Tab(hvg_layout(), label="Highly Variable Genes", tab_id="hvg_tab"),
						
						### CORRELATION NETWORK TAB ###
						dbc.Tab(correlation_tab_layout(), label="Correlations", tab_id="correlation_tab"),
						
					    ### ANNOTATION TAB ###
					    dbc.Tab(annotation_layout(), label="Exploration", tab_id="annotation_tab"),

						### DOWNLOAD ANALYSIS TAB ###
						dbc.Tab(exporting_layout(demo=demo), label="Save/export", tab_id="exporting_tab"),
					], id="main_tabs"),	
		    	], width=10),

		    	# Status panel
		    	dbc.Col(children=[
		    		status_layout()
		    	], width=2)
		    ])

		]),
	]) # end main layout
