import plotly.graph_objects as go
import numpy as np
from scipy import signal
import pandas as pd

######## C_in ########

# data = pd.read_csv('temp/cdata.csv')
# data.head()

# time = data['time'].values
# c_in = data['c_in'].values

# fig = go.Figure()

# fig.add_trace(go.Scatter(
#     x=time, y=c_in,
#     name='[Ca2+]_i',
#     marker_color='rgba(0, 0, 0, .7)' # gray
# ))

# fig.add_trace(go.Scatter(
#     x=QTDE_NCX, y=NCX02,
#     name='0.2 Hz',
#     marker_color='rgba(255, 102, 0, .9)' # Orange
# ))

# Set options common to all traces with fig.update_traces
# fig.update_traces(mode='lines', marker_line_width=2, marker_size=10)
# fig.update_layout(xaxis_title='<b>Time (s)<b>',
#                 yaxis_title='<b>Concentration (uM)<b>',
#                 yaxis_zeroline=False,
#                 xaxis_zeroline=False,
#                 font=dict(size=16),
#                 paper_bgcolor='rgba(0,0,0,0)',plot_bgcolor='rgba(0,0,0,0)')
# fig.update_yaxes(showline=True, linewidth=2, linecolor='black', showgrid=True, gridwidth=1, gridcolor='lightgray')
# fig.update_xaxes(showline=True, linewidth=2, linecolor='black', showgrid=True, gridwidth=1, gridcolor='lightgray')

# fig.show()

########## SNR ##########

# data = pd.read_csv('temp/SNR.csv')
# data.head()

# NCX01 = data['NCX 10'].values
# NCX02 = data['NCX 20'].values
# NCX03 = data['NCX 30'].values
# NCX04 = data['NCX 40'].values
# NCX06 = data['NCX 50'].values

# freq = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

# fig = go.Figure()

# fig.add_trace(go.Scatter(
#     x=freq, y=NCX01,
#     name='10 NCX',
#     marker_color='rgba(204, 0, 0, .9)' # Red
# ))

# fig.add_trace(go.Scatter(
#     x=freq, y=NCX02,
#     name='20 NCX',
#     marker_color='rgba(255, 102, 0, .9)' # Orange
# ))

# fig.add_trace(go.Scatter(
#     x=freq, y=NCX03,
#     name='30 NCX',
#     marker_color='rgba(255, 209, 26, .9)' # Yellow
# ))

# fig.add_trace(go.Scatter(
#     x=freq, y=NCX04,
#     name='40 NCX',
#     marker_color='rgba(0, 204, 0, .9)' # Green
# ))

# fig.add_trace(go.Scatter(
#     x=freq, y=NCX06,
#     name='50 NCX',
#     marker_color='rgba(0, 102, 255, .9)' # Blue
# ))

# # Set options common to all traces with fig.update_traces
# fig.update_traces(mode='markers+lines', marker_line_width=2, marker_size=10)
# fig.update_layout(xaxis_title='<b><i>f<i> (Hz)<b>',
#                 yaxis_title='<b>SNR (dB)<b>',
#                 yaxis_zeroline=False,
#                 xaxis_zeroline=False,
#                 font=dict(size=16),
#                 paper_bgcolor='rgba(0,0,0,0)',plot_bgcolor='rgba(0,0,0,0)')
# fig.update_yaxes(showline=True, linewidth=2, linecolor='black', showgrid=True, gridwidth=1, gridcolor='lightgray')
# fig.update_xaxes(showline=True, linewidth=2, linecolor='black', showgrid=True, gridwidth=1, gridcolor='lightgray')

# fig.show()


# t = np.linspace(0,200, 2000, endpoint=False)
# y = (signal.square(2*np.pi*0.6*t, duty=0.5)+1)/4

# fig = go.Figure()
# fig.add_trace(go.Scatter(x=t, y=y, line_color='rgb(0,0,0)'))

# fig.update_traces(hoverinfo='text+name', mode='lines')
# fig.update_layout(xaxis_title='t (s)', yaxis_title='[Ca]i (uM)', font=dict(size=14), paper_bgcolor='rgba(0,0,0,0)',plot_bgcolor='rgba(0,0,0,0)')
# fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='LightPink')
# # fig.update_xaxes(size=12)

# fig.show()

########## Semilog plot ##########

data = pd.read_csv('temp/cdata.csv')
data.head()

time_slot = data['Time Slot'].values
HL = data['Horizontal Line'].values
VL = data['Vertical Line'].values
C = data['Cross'].values
D = data['Diamond'].values

fig = go.Figure()

fig.add_trace(go.Scatter(
    x=time_slot, y=HL,
    name='Horizontal Line',
    marker_color='rgba(0, 102, 255, .9)' # Blue
))

fig.add_trace(go.Scatter(
    x=time_slot, y=VL,
    name='Vertical Line',
    marker_color='rgba(0, 204, 0, .9)' # Green
))

fig.add_trace(go.Scatter(
    x=time_slot, y=C,
    name='Cross',
    marker_color='rgba(0, 0, 0, .7)' # gray
))

fig.add_trace(go.Scatter(
    x=time_slot, y=D,
    name='Diamond',
    marker_color='rgba(204, 0, 0, .9)' # Red
))

fig.update_traces(mode='lines', marker_line_width=2, marker_size=10)
fig.update_layout(yaxis_type="log", xaxis_title='<b>Time Slot (s)<b>',
                yaxis_title='<b>BER<b>',
                yaxis_zeroline=False,
                xaxis_zeroline=False,
                font=dict(size=16),
                paper_bgcolor='rgba(0,0,0,0)',plot_bgcolor='rgba(0,0,0,0)')
fig.update_yaxes(showline=True, linewidth=2, linecolor='black', showgrid=True, gridwidth=1, gridcolor='lightgray')
fig.update_xaxes(showline=True, linewidth=2, linecolor='black', showgrid=True, gridwidth=1, gridcolor='lightgray')

fig.show()