
def integral(f,a,b,dt):
    num = a
    sum = 0
    while num <b:
        sum += f(num)*dt
        num = num + dt
    return sum

def EulerIntegral(f,a,b,n):
    if isinstance(n , int):
        if n <0:
            raise AttributeError
    else:
        raise AttributeError

    dt = (b-a)/n

    return integral(f,a,b,dt)