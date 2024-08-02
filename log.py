import os

class Log():
    def __init__(self, path="./log.txt"):
        self.path = path
        if os.path.exists(path):
            with open(path, 'w') as file:
                file.write("")
        else:
            open(path, 'w').close()

    def __call__(self, *args, seq=' ', end='\n'):
        file = open(self.path, 'a')
        for arg in args:
            arg = str(arg)
            print(arg, end=seq)
            file.write(arg)
            file.write(seq)
        print(end, end="")
        file.write(end)
        file.close()

    def changeDir(self, path="./log.txt"):
        self.path = path
        if os.path.exists(path):
            with open(path, 'w') as file:
                file.write("")
        else:
            open(path, 'w').close()

if __name__ == '__main__':
    a, b, c = 3.2, 4.5, 10086
    log = Log()
    log("gogo", "a:%.2f  b:%d"%(a, b), c)
    log("sdafsdf")
    log(a, b, c)

    log("code end~")