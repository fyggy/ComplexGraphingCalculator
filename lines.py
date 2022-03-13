from numpy import linspace, isnan, isinf, complex128, sign
from mpmath import fp

# TODO: verify if getters are nessessary
class Point:
    def __init__(self, input, output, derivative):
        self.input = input
        self.output = output
        self.derivative = derivative

    def __str__(self):
        return f"({self.input}, {self.output})"

class DeletedPoint(Point):
    input = None
    output = None
    derivative = None

class LinePart:
    def __init__(self, points):
        self.points = points

    @staticmethod
    def convert_single(z0, z1, d0, d1, direction):
        x0, y0 = z0.real, z0.imag
        x1, y1 = z1.real, z1.imag
        d0 *= direction
        d1 *= direction
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

    def convert(self, direction):
        if len(self.points) == 0:
            return [[]]
        output = [0] * (2 * len(self.points) - 1)
        next = self.points[0]
        z1 = next.output
        for i in range(len(self.points) - 1):
            current = self.points[i]
            next = self.points[i + 1]
            z0 = current.output
            z1 = next.output
            d0 = current.derivative
            d1 = next.derivative
            output[2 * i] = [z0.real, z0.imag]
            output[2 * i + 1] = list(LinePart.convert_single(z0, z1, d0, d1, direction))
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
            out = point.output
            dout = point.derivative
            if isinf(out) or isnan(out) or isinf(dout) or isnan(dout):
                self.points[i] = DeletedPoint

    @staticmethod
    def break_up(points):
        if len(points) == 0:
            return []

        output = []
        for i in range(len(points) - 1):
            current = points[i]
            next = points[i + 1]
            if current == DeletedPoint or next == DeletedPoint:
                start = i
                while current == DeletedPoint or next == DeletedPoint:
                    i += 1
                    current = points[i]
                    try:
                        next = points[i + 1]
                    except IndexError:
                        break
                if start == 0:
                    output.append(LinePart([]))
                else:
                    output.append(LinePart(points[:start+1]))
                output += Line.break_up(points[i+1:])
                return output

            next = points[i + 1]
            diff = next.input - current.input
            s = diff


            check = ((abs(next.output - (current.output + (s * current.derivative)))) ** 2) / (abs(current.output) + s)
            if abs(check) >= (abs(current.derivative)):
                print("broken")
                output.append(LinePart(points[:i+1]))
                output += Line.break_up(points[i+1:])
                return output

        output.append(LinePart(points))
        return output

    def break_fully(self):
        self.lineparts = Line.break_up(self.points)

    def convert(self, direction):
        output = []
        for i in self.lineparts:
            output.append(i.convert(direction))

        return output
