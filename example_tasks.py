import numpy as np
import plotly as py
import plotly.graph_objects as go


class ExampleTasks:

    def task1(self, a=0, b=0, u_0=0, n=0):

        if not all([u_0 > 0, n > 0]):
            return

        delta_t = (b - a) / (n - 1)

        u = np.array([u_0])
        x = np.array([0])

        for i in range(0, n - 1):
            u = np.append(u, round(u[i] * (1 + delta_t) + 0.1 * np.exp(x[i]) * delta_t, 2))
            x = np.append(x, round(x[i] + delta_t, 2))

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=x, y=u, mode='lines+markers', name='Numerische Lösung'))
        fig.add_trace(go.Scatter(x=x, y=(1 + 0.1 * x) * np.exp(x), mode='lines+markers', name='Exakte Lösung'))
        fig.update_xaxes(title_text='x')
        fig.update_yaxes(title_text='y')
        fig.update_layout(title='EulerVorwaerts 1')
        py.offline.plot(fig, config={'scrollZoom': True}, filename='Task1.html')

    def task2(self, a=0, b=0, u_0=0, n=0):

        if not all([u_0 > 0, n > 0]):
            return

        delta_t = (b - a) / (n - 1)

        u = [u_0]
        x = np.linspace(a, b, n)

        for i in range(0, n - 1):
            u = np.append(u, round(u[i] * (1 + delta_t) + 0.1 * np.exp(x[i]) * delta_t, 2))

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=x, y=u, mode='lines+markers', name='Numerische Lösung'))
        fig.add_trace(go.Scatter(x=x, y=(1 + 0.1 * x) * np.exp(x), mode='lines+markers', name='Exakte Lösung'))
        fig.update_xaxes(title_text='x')
        fig.update_yaxes(title_text='y')
        fig.update_layout(title='EulerVorwaerts 2')
        fig.write_image('Task2.pdf')

    def task3(self, l=0, t=0, c=0, i_max=0, delta_t=0):
        fig_0 = None
        fig_3 = go.Figure()

        if not all([l > 0, t > 0, c > 0, i_max > 0, delta_t > 0]):
            return

        x = np.linspace(0, l, i_max)
        delta_x = l / (i_max - 1)
        n_max = int(t / delta_t) + 1

        u = np.zeros((i_max, n_max))
        v = np.zeros((i_max, n_max))
        w = np.zeros((n_max, n_max))

        for i in range(0, n_max):
            w[i, 0] = i * delta_t

        for i in range(0, i_max):
            u[i, 0] = round(np.exp(-1 * np.square(x[i] - 4)), 4)
            v[i, 0] = u[i, 0]

        for i in range(1, n_max):
            u[0, i] = u[i_max - 1, i - 1]

            for j in range(1, i_max):
                u[j, i] = round(u[j, i - 1] - delta_t / delta_x * c * (u[j, i - 1] - u[j - 1, i - 1]), 4)

            tmpl = round(c * w[i, 0] + 4, 4)

            for k in range(0, i_max):
                v[k, i] = round(np.exp(-1 * np.square(x[k] - tmpl)), 4)

        for i in range(0, n_max):
            fig_0 = go.Figure(
                data=[go.Scatter(x=x, y=u[0:i_max, 0])],
                frames=[go.Frame(data=[go.Scatter(x=x, y=u[0:i_max, i])])])
        fig_0.update_xaxes(title_text='x')
        fig_0.update_yaxes(title_text='y')
        py.offline.plot(fig_0, config={'scrollZoom': True}, filename='Task3_plot0.html')

        fig_1 = go.Figure(data=[go.Surface(z=u)])
        py.offline.plot(fig_1, config={'scrollZoom': True}, filename='Task3_plot1.html')
        fig_1.write_image('Task3_plot1.pdf')

        fig_2 = go.Figure(data=[go.Contour(z=u)])
        py.offline.plot(fig_2, config={'scrollZoom': True}, filename='Task3_plot2.html')

        for i in range(0, n_max):
            fig_3.add_trace(go.Scatter(x=x, y=u[0:i_max, i], line=dict(color='firebrick')))
            fig_3.add_trace(go.Scatter(x=x + 0.1, y=v[0:i_max, i], line=dict(color='royalblue')))
        fig_3.update_xaxes(title_text='x')
        fig_3.update_yaxes(title_text='y')
        fig_3.update_layout(showlegend=False)
        py.offline.plot(fig_3, config={'scrollZoom': True}, filename='Task3_plot3.html')

    def task4(self, l=0, t=0, c=0, i_max=0, cfl=0):

        fig_0 = None
        fig_3 = go.Figure()

        if not all([l > 0, t > 0, c != 0, i_max > 0, cfl > 0]):
            return

        x = np.linspace(0, l, i_max)
        delta_x = l / (i_max - 1)
        delta_t = cfl * delta_x / np.abs(c)
        n_max = int(np.ceil(t / delta_t)) + 1

        c_plus = np.amax(c, initial=0)
        c_minus = np.amin(c, initial=0)

        u = np.zeros((i_max, n_max))
        v = np.zeros((i_max, n_max))
        w = np.zeros((n_max, n_max))  # t bei Prof. Dr. Voß

        for i in range(0, n_max):
            w[i, 0] = i * delta_t

        for i in range(0, i_max):
            u[i, 0] = round(np.exp(-1 * np.square(x[i] - 4)), 4)
            v[i, 0] = u[i, 0]

        for i in range(1, n_max):
            if c > 0:
                u[0, i] = u[i_max - 1, i - 1]
                tmp_plus = delta_t / delta_x * c_plus * (u[i_max - 1, i - 1] - u[i_max - 2, i - 1])
                u[i_max - 1, i] = round(u[i_max - 1, i - 1] - tmp_plus, 4)

            if c < 0:
                u[i_max - 1, i] = u[0, i - 1]
                tmp_minus = delta_t / delta_x * c_minus * (u[1, i - 1] - u[0, i - 1])
                u[0, i] = round(u[0, i - 1] - tmp_minus, 4)

            for j in range(1, i_max - 1):
                tmp_plus = delta_t / delta_x * c_plus * (u[j, i - 1] - u[j - 1, i - 1])
                tmp_minus = delta_t / delta_x * c_minus * (u[j + 1, i - 1] - u[j, i - 1])
                u[j, i] = round(u[j, i - 1] - tmp_plus - tmp_minus, 4)

            tmpl = round(c * w[i, 0] + 4, 4)

            for k in range(0, i_max):
                v[k, i] = round(np.exp(-1 * np.square(x[k] - tmpl)), 4)

        for i in range(0, n_max):
            fig_0 = go.Figure(
                data=[go.Scatter(x=x, y=u[0:i_max, 0])],
                frames=[go.Frame(data=[go.Scatter(x=x, y=u[0:i_max, i])])])
        fig_0.update_xaxes(title_text='x')
        fig_0.update_yaxes(title_text='y')
        py.offline.plot(fig_0, config={'scrollZoom': True}, filename='Task4_plot0.html')

        fig_1 = go.Figure(data=[go.Surface(z=u)])
        py.offline.plot(fig_1, config={'scrollZoom': True}, filename='Task4_plot1.html')
        fig_1.write_image('Task4_plot1.pdf')

        fig_2 = go.Figure(data=[go.Contour(z=u)])
        py.offline.plot(fig_2, config={'scrollZoom': True}, filename='Task4_plot2.html')

        for i in range(0, n_max):
            fig_3.add_trace(go.Scatter(x=x, y=u[0:i_max, i], line=dict(color='firebrick')))
            fig_3.add_trace(go.Scatter(x=x + 0.1, y=v[0:i_max, i], line=dict(color='royalblue')))
        fig_3.update_xaxes(title_text='x')
        fig_3.update_yaxes(title_text='y')
        fig_3.update_layout(showlegend=False)
        py.offline.plot(fig_3, config={'scrollZoom': True}, filename='Task4_plot3.html')

    def task5(self, l=0, t=0, c=0, i_max=0, cfl=0):

        if not all([l > 0, t > 0, c != 0, i_max > 0, cfl > 0]):
            return

        fig_0 = None
        fig_3 = go.Figure()

        x = np.linspace(0, l, i_max)
        delta_x = l / (i_max - 1)
        delta_t = cfl * delta_x / np.abs(c)
        n_max = int(np.ceil(t / delta_t)) + 1

        u = np.zeros((i_max, n_max))
        v = np.zeros((i_max, n_max))
        w = np.linspace(0, t, n_max)  # t bei Prof. Dr. Voß

        for i in range(0, i_max):
            if x[i] < 2:
                u[i, 0] = 1
            elif x[i] < 3:
                u[i, 0] = round(0.15 * np.sin(2 * np.pi * (x[i] - 2)) + 1, 4)
            else:
                u[i, 0] = 1

        for n in range(1, n_max):
            u[0, n] = 1
            u[i_max - 1, n] = 1

            for i in range(1, i_max - 1):
                c_plus = np.amax(c, initial=0)
                c_minus = np.amin(c, initial=0)
                tmp_plus = round(delta_t / delta_x * c_plus * (u[i, n - 1] - u[i - 1, n - 1]), 4)
                tmp_minus = round(delta_t / delta_x * c_minus * (u[i + 1, n - 1] - u[i, n - 1]), 4)
                u[i, n] = round(u[i, n - 1] - tmp_plus - tmp_minus, 4)

        for i in range(0, n_max):
            fig_0 = go.Figure(
                data=[go.Scatter(x=x, y=u[0:i_max, 0])],
                frames=[go.Frame(data=[go.Scatter(x=x, y=u[0:i_max, i])])])
        fig_0.update_xaxes(title_text='x')
        fig_0.update_yaxes(title_text='y')
        py.offline.plot(fig_0, config={'scrollZoom': True}, filename='Task4_plot0.html')

        fig_1 = go.Figure(data=[go.Surface(z=u)])
        py.offline.plot(fig_1, config={'scrollZoom': True}, filename='Task4_plot1.html')
        fig_1.write_image('Task4_plot1.pdf')

        fig_2 = go.Figure(data=[go.Contour(z=u)])
        py.offline.plot(fig_2, config={'scrollZoom': True}, filename='Task4_plot2.html')

        for i in range(0, n_max):
            fig_3.add_trace(go.Scatter(x=x, y=u[0:i_max, i], line=dict(color='firebrick')))
            fig_3.add_trace(go.Scatter(x=x + 0.1, y=v[0:i_max, i], line=dict(color='royalblue')))
        fig_3.update_xaxes(title_text='x')
        fig_3.update_yaxes(title_text='y')
        fig_3.update_layout(showlegend=False)
        py.offline.plot(fig_3, config={'scrollZoom': True}, filename='Task4_plot3.html')
