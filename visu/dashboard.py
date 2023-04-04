import os
import pandas as pd
import plotly.express as px


path_data = ''

sub = []
ses = []
run = []

df = pd.read_csv(os.path.join(path_data, sub, ses, f'{run}_data.csv'))

signal_type = []

fig = px.line(df, x='time', y=signal_type, title=f"{signal_type} timeserie : {sub} {ses} {run}")

fig.update_xaxes(
    rangeslider_visible=True,
    rangeselector=dict(
        buttons=list([
            dict(step="all")
        ])
    )
)

fig.show()