# https://github.com/jupyter-widgets/ipywidgets/issues/2522
c.NotebookApp.tornado_settings = {"websocket_max_message_size": 100 * 1024 * 1024}

