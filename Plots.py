#plots.py>

# plotagens orbitais #
def plot_animacao_orbita(dataframe, size, numero_divisoes):
    import pandas as pd
    import plotly.express as px
    from icosphere import icosphere
    import plotly.graph_objects as go
    df = dataframe
    size = size
    fig = px.scatter_3d(df, x="X_ECI", y="Y_ECI", z="Z_ECI",
                        range_x=[-size, size],
                        range_y=[-size, size],
                        range_z=[-size, size],
                        animation_group="Data",
                        animation_frame="Tempo",
                        title="Spaceflight",
                        size_max=100, opacity=0.7)


    fig.show()
    return

def plot_groundtrack_3D(dataframe):
    import plotly.graph_objects as go
    import pandas as pd

    df = dataframe

    scl = ['rgb(213,62,79)', 'rgb(244,109,67)', 'rgb(253,174,97)', \
           'rgb(254,224,139)', 'rgb(255,255,191)', 'rgb(230,245,152)', \
           'rgb(171,221,164)', 'rgb(102,194,165)', 'rgb(50,136,189)'
           ]
    n_colors = len(scl)

    fig = go.Figure()


    fig.add_trace(go.Scattergeo(
        lon=df['longitude'],
        lat=df['latitude'],
        mode='lines',
        line=dict(width=2, color='rgb(213,62,79)'
                      ),
        connectgaps=False)
        )

    fig.update_layout(
        title_text='Contour lines over globe<br>(Click and drag to rotate)',
        showlegend=False,
        geo=dict(
            showland=True,
            showcountries=True,
            showocean=True,
            countrywidth=0.5,
            landcolor='rgb(34,139,34)',
            lakecolor='rgb(0, 255, 255)',
            oceancolor='rgb(0,0,139)',
            projection=dict(
                type='orthographic',
                rotation=dict(
                    lon=-100,
                    lat=40,
                    roll=0
                ),

            ),
            lonaxis=dict(
                showgrid=True,
                gridcolor='rgb(0,0,0)',
                gridwidth=0.5
            ),
            lataxis=dict(
                showgrid=True,
                gridcolor='rgb(0,0,0)',
                gridwidth=0.5
            )
        )
    )

    return fig.show()

def plot_groundtrack_2D(dataframe):


    import plotly.graph_objects as go
    import plotly.io as pio
    from PIL import Image
    df = dataframe
    fig = go.Figure(data=go.Scattergeo(
        lat=df['latitude'],
        lon=df['longitude'],
        mode='lines',
        line=dict(width=2, color='red'),
        name='A1',
        showlegend=True
    ))

    fig.update_layout(
        title= dict(font_color = 'red',
                    text = 'Groundtrack 2D',
                    x = 0.5,
                    ),
        showlegend=True,
        geo=dict(
            showland=True,
            showcountries=True,
            showocean=True,
            countrywidth=0.5,
            landcolor='rgb(255,255,255)',
            lakecolor='rgb(240,248,255)',
            oceancolor='rgb(240,248,255)',
            projection_type="equirectangular",
            coastlinewidth=1,
            lataxis=dict(
                dtick=30,
                gridcolor='rgb(0, 0, 0)',
                griddash = "dash",
                gridwidth=0.5,
                range=[-90, 90],
                showgrid=True,
                tick0 = -90

            ),
            lonaxis=dict(
                range=[-180, 180],
                showgrid=True,
                gridcolor='rgb(0, 0, 0)',
                griddash="dash",
                gridwidth=0.5,
                dtick=60,
                tick0 = -180
            )

        )
    )

    fig.show()
    return
def plot_ground2d_novo(dataframe):
    import plotly.graph_objects as go
    import plotly.io as pio
    from PIL import Image

    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    df = dataframe
    f = plt.figure(figsize=(10, 7.5))
    m = Basemap(projection="mill", lon_0=0);
    m.drawcoastlines()
    m.drawparallels(np.arange(-90, 91, 30), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180, 181, 60), labels=[0, 0, 0, 1])

    plt.savefig('groundtrack.png')
    # Create figure
    fig = go.Figure()

    pyLogo = Image.open("groundtrack.png")

    fig = go.Figure()

    # Constants
    img_width = 1600
    img_height = 900
    scale_factor = 0.5

    # Add invisible scatter trace.
    # This trace is added to help the autoresize logic work.


    # Configure axes
    fig.update_xaxes(
        visible=False,
        range=[0, img_width * scale_factor]
    )

    fig.update_yaxes(
        visible=False,
        range=[0, img_height * scale_factor],
        # the scaleanchor attribute ensures that the aspect ratio stays constant
        scaleanchor="x"
    )

    # Add image
    fig.add_layout_image(
        dict(
            x=0,
            sizex=img_width * scale_factor,
            y=img_height * scale_factor,
            sizey=img_height * scale_factor,
            xref="x",
            yref="y",
            opacity=1.0,
            layer="below",
            sizing="stretch",
            source=pyLogo)
    )

    # Configure other layout
    fig.update_layout(
        width=img_width * scale_factor,
        height=img_height * scale_factor,
        margin={"l": 0, "r": 0, "t": 0, "b": 0},
    )
    fig.add_trace(
        go.Scatter(
            x=df['longitude'],
            y=df['latitude'],
            mode="markers",
            marker_opacity=0
        )
    )
    # Disable the autosize on double click because it adds unwanted margins around the image
    # More detail: https://plotly.com/python/configuration-options/
    fig.show(config={'doubleClick': 'reset'})
# Plotagem termica #

def calor_solar(dataframe):
    import plotly.express as px
    Q_sol = dataframe
    linhas = px.line(Q_sol, y=['Solar 1', 'Solar 2', 'Solar 3', 'Solar 4', 'Solar 5', 'Solar 6'])

    linhas.update_layout(title=dict(font=dict(color = 'black'),
                                    text = 'Radiação Solar',
                                    x = 0.5),
                         showlegend = True,
                         legend = dict(bgcolor = 'rgb(228, 232, 255)',
                                       bordercolor = 'rgb(0, 0, 0)',
                                       borderwidth = 1.0,
                                       itemdoubleclick = "toggleothers",
                                       title_text = 'Face'
                                       ),
                         xaxis=dict(
                             showgrid=True,
                             gridcolor='rgb(255, 255, 255)',
                             gridwidth=1,
                             dtick=len(Q_sol)*0.1
                         ),
                         yaxis=dict(
                             showgrid=True,
                             gridcolor='rgb(255, 255, 255)',
                             gridwidth=1,
                             dtick=100
                         ),
                         autosize=True,
                         plot_bgcolor = 'rgb(228, 232, 255)',
                      xaxis_title='Posição Orbital',
                      yaxis_title='Radiação [W/m^2]')
    linhas.show()

    return

def calor_albedo(dataframe):
    import plotly.express as px
    Q_alb = dataframe
    linhas = px.line(Q_alb, y=['Albedo 1', 'Albedo 2', 'Albedo 3', 'Albedo 4', 'Albedo 5', 'Albedo 6'])
    linhas.update_layout(title=dict(font=dict(color = 'black'),
                                    text = 'Radiação de Albedo',
                                    x = 0.5),
                         showlegend = True,
                         legend = dict(bgcolor = 'rgb(228, 232, 255)',
                                       bordercolor = 'rgb(0, 0, 0)',
                                       borderwidth = 1.0,
                                       itemdoubleclick = "toggleothers",
                                       title_text = 'Face'
                                       ),
                         xaxis=dict(
                             showgrid=True,
                             gridcolor='rgb(255, 255, 255)',
                             gridwidth=1,
                             dtick=100
                         ),
                         yaxis=dict(
                             showgrid=True,
                             gridcolor='rgb(255, 255, 255)',
                             gridwidth=1,
                             dtick=100
                         ),
                         autosize=True,
                         plot_bgcolor = 'rgb(228, 232, 255)',
                      xaxis_title='Posição Orbital',
                      yaxis_title='Radiação [W/m^2]')
    linhas.show()
    return

def calor_IR_Terra(dataframe):

    import plotly.express as px
    Q_ir = dataframe
    linhas = px.line(Q_ir, y = ['IR Terra 1', 'IR Terra 2', 'IR Terra 3', 'IR Terra 4', 'IR Terra 5', 'IR Terra 6'])
    linhas.update_layout(title=dict(font=dict(color = 'black'),
                                    text = 'Radiação IR Terra',
                                    x = 0.5),
                         showlegend = True,
                         legend = dict(bgcolor = 'rgb(228, 232, 255)',
                                       bordercolor = 'rgb(0, 0, 0)',
                                       borderwidth = 1.0,
                                       itemdoubleclick = "toggleothers",
                                       title_text = 'Face'
                                       ),
                         xaxis=dict(
                             showgrid=True,
                             gridcolor='rgb(255, 255, 255)',
                             gridwidth=1,
                             dtick=100
                         ),
                         yaxis=dict(
                             showgrid=True,
                             gridcolor='rgb(255, 255, 255)',
                             gridwidth=1,
                             dtick=100
                         ),
                         autosize=True,
                         plot_bgcolor = 'rgb(228, 232, 255)',
                      xaxis_title='Posição Orbital',
                      yaxis_title='Radiação [W/m^2]')
    linhas.show()
    return

def calor_total(dataframe):

    import plotly.express as px
    Q_total = dataframe
    linhas = px.line(Q_total, y=['Total 1', 'Total 2', 'Total 3', 'Total 4', 'Total 5', 'Total 6'])
    linhas.update_layout(title=dict(font=dict(color = 'black'),
                                    text = 'Radiação Total em Cada Face',
                                    x = 0.5),
                         showlegend = True,
                         legend = dict(bgcolor = 'rgb(228, 232, 255)',
                                       bordercolor = 'rgb(0, 0, 0)',
                                       borderwidth = 1.0,
                                       itemdoubleclick = "toggleothers",
                                       title_text = 'Face'
                                       ),
                         xaxis=dict(
                             showgrid=True,
                             gridcolor='rgb(255, 255, 255)',
                             gridwidth=1,
                             dtick=100
                         ),
                         yaxis=dict(
                             showgrid=True,
                             gridcolor='rgb(255, 255, 255)',
                             gridwidth=1,
                             dtick=100
                         ),
                         autosize=True,
                         plot_bgcolor = 'rgb(228, 232, 255)',
                      xaxis_title='Posição Orbital',
                      yaxis_title='Radiação [W/m^2]')
    linhas.show()

    return
