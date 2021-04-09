from example_tasks import ExampleTasks


def main():
    tasks = ExampleTasks()
    tasks.task1(b=1, u_0=1, n=11)
    tasks.task2(b=1, u_0=1, n=11)
    tasks.task3(l=16, t=10, c=1, i_max=17, delta_t=1)
    tasks.task4(l=16, t=10, c=1, i_max=100, cfl=0.9)
    tasks.task5(l=4, t=1, c=1, i_max=400, cfl=0.9)


if __name__ == '__main__':
    main()