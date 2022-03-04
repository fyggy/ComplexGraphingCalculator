from numpy import linspace, isnan, isinf, complex128

# TODO: verify if getters are nessessary
class Point:
    def __init__(self, input, output, derivative):
        self.input = input
        self.output = output
        self.derivative = derivative

    def get_input(self):
        return self.input

    def get_output(self):
        return self.output

    def get_derivative(self):
        return self.derivative

    def __str__(self):
        return f"({self.input}, {self.output})"
        
class LinePart:
    def __init__(self, points):
        self.points = points

    @staticmethod
    def convert_single(z0, z1, d0, d1):
        x0, y0 = z0.real, z0.imag
        x1, y1 = z1.real, z1.imag

        if d0.real == 0 and d1.real == 0:
            tmp = (z0 + z1) / 2
            return (tmp.real, tmp.imag)

        elif d0.real == 0:
            m1 = d1.imag / d1.real
            return (x0, m1 * (x0 - x1) + y1)

        elif d1.real == 0:
            m0 = d0.imag / d0.real
            return (x1, m0 * (x1 - x0) + y0)

        else:
            m0 = d0.imag / d0.real
            m1 = d1.imag / d1.real

            if m0 == m1:
                tmp = (z0 + z1) / 2
                return (tmp.real, tmp.imag)

            else:
                tmp = (m0 * x0 - m1 * x1 + y1 - y0) / (m0 - m1)
                return (tmp, m0 * (tmp - x0) + y0)

    def convert(self):
        output = [0] * (2 * len(self.points) + 1)
        next = self.points[0]
        for i in range(len(self.points) - 1):
            current = self.points[i]
            next = self.points[i + 1]
            z0 = current.get_output()
            z1 = next.get_output()
            d0 = current.get_derivative()
            d1 = next.get_derivative()
            output[2 * i] = [z0.real, z0.imag]
            output[2 * i + 1] = list(LinePart.convert_single(z0, z1, d0, d1))

        output[-1] = [z1.real, z1.imag]
        return output

class Line:
    def __init__(self, start, end, step, function, dfunction):
        self.start = start
        self.end = end
        self.step = step
        self.function = function
        self.dfunction = dfunction
        self.input = linspace(start, end, num=step, dtype=complex128)

    def calculate(self):
        output = self.function(self.input)
        doutput = self.dfunction(self.input)
        self.points = [Point] * len(self.input)
        for i, (inp, out, dout) in enumerate(zip(self.input, output, doutput)):
            self.points[i] = Point(inp, out, dout)

    def trim(self):
        for i, point in enumerate(self.points):
            out = point.get_output()
            dout = point.get_derivative()
            if isinf(out) or isnan(out) or isinf(dout) or isnan(dout):
                del self.points[i]

    @staticmethod
    def break_up(points):
        output = []
        for i in range(len(points) - 1):
            current = points[i]
            next = points[i + 1]
            diff = next.get_input() - current.get_input()

            #TODO: make this work properly
            if abs(diff) > 0.1:
                output.append(LinePart(points[:i]))
                output += Line.break_up(points[i+1:])
                return output

            else:
                s = abs(diff)
                check = ((abs(next.get_output() - current.get_output() - (s * current.get_derivative()))) ** 2) / (abs(current.get_output()) + s)
                if check >= 0.1:
                    output.append(LinePart(points[:i]))
                    print(points[i+1:])
                    output += Line.break_up(points[i+1:])
                    print(points)
                    return output
                else:
                    output.append(LinePart(points))
                    return output

    def break_fully(self):
        self.lineparts = Line.break_up(self.points)

    def convert(self):
        output = []
        for i in self.lineparts:
            output.append(i.convert())

        return output
