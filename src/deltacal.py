__author__ = 'Gian Paolo Jesi'

from argparse import ArgumentParser
from math import cos, sin, pi, sqrt, isnan
from numpy import zeros, empty, array

SUPPORTED_FIRMWARES = ['Marlin', 'MarlinRC', 'RRFW', 'Smoothieware', 'Repetier']
firmware = 'RRFW'
degreesToRadians = pi / 180.0


class DeltaParameters(object):
    def __init__(self, diagonal=215.0, delta_radius=105.6, height=250.0, xstop=0.0, ystop=0.0, zstop=0.0, xadj=0.0,
                 yadj=0.0, zadj=0.0):
        self.diagonal = diagonal
        self.delta_radius = delta_radius
        self.homedHeight = height
        self.xstop = xstop
        self.ystop = ystop
        self.zstop = zstop
        self.xadj = xadj
        self.yadj = yadj
        self.zadj = zadj

        self.Xbc = 0
        self.Xca = 0
        self.Xab = 0
        self.Ybc = 0
        self.Yca = 0
        self.Yab = 0

        self.coreFa = 0
        self.coreFb = 0
        self.coreFc = 0

        self.Q = 0
        self.Q2 = 0
        self.D2 = 0

        self.towerX = empty(3)
        self.towerY = empty(3)

        self.homedCarriageHeight = 0

        self.recalc()

    def transform(self, machinePos, axis):
        return machinePos[2] + sqrt(
            self.D2 - ((machinePos[0] - self.towerX[axis]) * (machinePos[0] - self.towerX[axis])) -
            ((machinePos[1] - self.towerY[axis]) * (machinePos[1] - self.towerY[axis])))

    def inverse_transform(self, Ha, Hb, Hc):
        Fa = self.coreFa + (Ha * Ha)
        Fb = self.coreFb + (Hb * Hb)
        Fc = self.coreFc + (Hc * Hc)

        # Setup PQRSU such that x = -(S - uz)/P, y = (P - Rz)/Q
        P = (self.Xbc * Fa) + (self.Xca * Fb) + (self.Xab * Fc)
        S = (self.Ybc * Fa) + (self.Yca * Fb) + (self.Yab * Fc)

        R = 2 * ((self.Xbc * Ha) + (self.Xca * Hb) + (self.Xab * Hc))
        U = 2 * ((self.Ybc * Ha) + (self.Yca * Hb) + (self.Yab * Hc))

        R2 = R * R
        U2 = U * U

        A = U2 + R2 + self.Q2
        minus_half_b = S * U + P * R + Ha * self.Q2 + self.towerX[0] * U * self.Q - self.towerY[0] * R * self.Q
        C = ((S + self.towerX[0] * self.Q) * (S + self.towerX[0] * self.Q)) + \
            ((P - self.towerY[0] * self.Q) * (P - self.towerY[0] * self.Q)) + \
            ((Ha * Ha) - self.D2) * self.Q2

        rslt = (minus_half_b - sqrt((minus_half_b * minus_half_b) - A * C)) / A
        if isnan(rslt):
            raise Exception(
                "At least one probe point is not reachable. Please correct your delta radius, diagonal rod length, or probe coordinates.")

        return rslt

    def recalc(self):
        self.towerX[0] = -(self.delta_radius * cos((30 + self.xadj) * degreesToRadians))
        self.towerY[0] = -(self.delta_radius * sin((30 + self.xadj) * degreesToRadians))
        self.towerX[1] = +(self.delta_radius * cos((30 - self.yadj) * degreesToRadians))
        self.towerY[1] = -(self.delta_radius * sin((30 - self.yadj) * degreesToRadians))
        self.towerX[2] = -(self.delta_radius * sin(self.zadj * degreesToRadians))
        self.towerY[2] = +(self.delta_radius * cos(self.zadj * degreesToRadians))

        self.Xbc = self.towerX[2] - self.towerX[1]
        self.Xca = self.towerX[0] - self.towerX[2]
        self.Xab = self.towerX[1] - self.towerX[0]
        self.Ybc = self.towerY[2] - self.towerY[1]
        self.Yca = self.towerY[0] - self.towerY[2]
        self.Yab = self.towerY[1] - self.towerY[0]
        self.coreFa = (self.towerX[0] * self.towerX[0]) + (self.towerY[0] * self.towerY[0])
        self.coreFb = (self.towerX[1] * self.towerX[1]) + (self.towerY[1] * self.towerY[1])
        self.coreFc = (self.towerX[2] * self.towerY[2]) + (self.towerY[2] * self.towerY[2])
        self.Q = 2 * (self.Xca * self.Yab - self.Xab * self.Yca)
        self.Q2 = self.Q * self.Q
        self.D2 = self.diagonal * self.diagonal

        # Calculate the base carriage height when the printer is homed.
        temp_height = self.diagonal  # any sensible height will do here, probably even zero
        self.homedCarriageHeight = self.homedHeight + temp_height - self.inverse_transform(temp_height, temp_height,
                                                                                           temp_height)

    def computeDerivative(self, deriv, ha, hb, hc):
        perturb = 0.2  # perturbation amount in mm or degrees
        hiParams = DeltaParameters(self.diagonal, self.delta_radius, self.homedHeight, self.xstop, self.ystop,
                                   self.zstop, self.xadj, self.yadj, self.zadj)
        loParams = DeltaParameters(self.diagonal, self.delta_radius, self.homedHeight, self.xstop, self.ystop,
                                   self.zstop, self.xadj, self.yadj, self.zadj)

        if deriv == 0 or deriv == 1 or deriv == 2:
            pass

        elif deriv == 3:
            hiParams.delta_radius += perturb
            loParams.delta_radius -= perturb

        elif deriv == 4:
            hiParams.xadj += perturb
            loParams.xadj -= perturb

        elif deriv == 5:
            hiParams.yadj += perturb
            loParams.yadj -= perturb

        elif deriv == 6:
            hiParams.diagonal += perturb
            loParams.diagonal -= perturb

        else:
            print "Error!"

        hiParams.recalc()
        loParams.recalc()

        ha_hi = ha + perturb if (deriv == 0) else ha
        hb_hi = hb + perturb if (deriv == 1) else hb
        hc_hi = hc + perturb if (deriv == 2) else hc

        ha_lo = ha - perturb if (deriv == 0) else ha
        hb_lo = hb - perturb if (deriv == 1) else hb
        hc_lo = hc - perturb if (deriv == 2) else hc

        zHi = hiParams.inverse_transform(ha_hi, hb_hi, hc_hi)
        zLo = loParams.inverse_transform(ha_lo, hb_lo, hc_lo)

        return (zHi - zLo) / (2 * perturb)

    def normalise_end_stop_adjustments(self):
        eav = min(self.xstop, min(self.ystop, self.zstop)) if (
            firmware == "Marlin" or firmware == "MarlinRC" or firmware == "Repetier") else (
                                                                                               self.xstop + self.ystop + self.zstop) / 3.0
        self.xstop -= eav
        self.ystop -= eav
        self.zstop -= eav

        self.homedHeight += eav
        self.homedCarriageHeight += eav  # no need for a full recalc, self is sufficient

    def adjust(self, num_factors, v, norm=True):
        old_carriage_height = self.homedCarriageHeight + self.xstop  # save for later  # Update endstop adjustments
        self.xstop += v[0]
        self.ystop += v[1]
        self.zstop += v[2]

        if norm:
            self.normalise_end_stop_adjustments()

        if num_factors >= 4:
            self.delta_radius += v[3]

            if num_factors >= 6:
                self.xadj += v[4]
                self.yadj += v[5]

                if num_factors == 7:
                    self.diagonal += v[6]

            self.recalc()

        # Adjusting the diagonal and the tower positions affects the homed carriage height.
        # We need to adjust homedHeight to allow for self, to get the change that was requested in the endstop corrections.
        height_error = self.homedCarriageHeight + self.xstop - old_carriage_height - v[0]
        self.homedHeight -= height_error
        self.homedCarriageHeight -= height_error


def swap_rows(A, i, j, num_cols=-1):
    if i != j:
        for k in range(0, num_cols):
            temp = A[i][k]
            A[i][k] = A[j][k]
            A[j][k] = temp


def gaussjordan(A, solution, num_rows):
    """Gauss Jordan elimination algorithm. Straight conversion from DC42 javascript code.

    :param A: matrix
    :param solution: solution array. It is filled by the function
    :param num_rows: number of rows
    :return: None
    """
    for i in range(0, num_rows):
        # Swap the rows around for stable Gauss-Jordan elimination
        vmax = abs(A[i][i])
        for j in range(i + 1, num_rows):
            rmax = abs(A[j][i])
            if rmax > vmax:
                print "swapping"
                swap_rows(A, i, j)
                vmax = rmax

        # Use row i to eliminate the ith element from previous and subsequent rows
        v = A[i][i]
        for j in range(0, i):
            factor = A[j][i] / v
            A[j][i] = 0.0
            for k in range(i + 1, num_rows + 1):
                A[j][k] -= A[i][k] * factor

        for j in range(i + 1, num_rows):
            factor = A[j][i] / v
            A[j][i] = 0.0
            for k in range(i + 1, num_rows + 1):
                A[j][k] -= A[i][k] * factor

    for i in range(0, num_rows):
        solution[i] = (A[i][num_rows] / A[i][i])


def do_delta_calibration(deltapar, numFactors=6, numPoints=10, xBedProbePoints=[], yBedProbePoints=[],
                         zBedProbePoints=[], normalise=True):
    if numFactors != 3 and numFactors != 4 and numFactors != 6 and numFactors != 7:
        raise Exception("Error: " + numFactors + " factors requested but only 3, 4, 6 and 7 supported")

    if numFactors > numPoints:
        raise Exception("Error: need at least as many points as factors you want to calibrate")

    # Transform the probing points to motor endpoints and store them in a matrix,
    # so that we can do multiple iterations using the same data
    probe_motor_positions = empty(shape=(numPoints, 3))
    corrections = zeros(numPoints)
    initial_sum_of_squares = 0.0
    for i in range(0, numPoints):
        xp = xBedProbePoints[i]
        yp = yBedProbePoints[i]
        machine_pos = array([xp, yp, 0.0])
        probe_motor_positions[i][0] = deltapar.transform(machine_pos, 0)
        probe_motor_positions[i][1] = deltapar.transform(machine_pos, 1)
        probe_motor_positions[i][2] = deltapar.transform(machine_pos, 2)

        initial_sum_of_squares += zBedProbePoints[i] * zBedProbePoints[i]

    # print "Motor positions: ", probe_motor_positions

    # Do 1 or more Newton-Raphson iterations
    iteration = 0
    while (True):
        # Build a Nx7 matrix of derivatives with respect to xa, xb, yc, za, zb, zc, diagonal.
        # derivative_matrix = Matrix(numPoints, numFactors)
        derivative_matrix = empty(shape=(numPoints, numFactors))
        for i in range(numPoints):
            for j in range(numFactors):
                derivative_matrix[i][j] = deltapar.computeDerivative(j, probe_motor_positions[i][0],
                                                                     probe_motor_positions[i][1],
                                                                     probe_motor_positions[i][2])

        # print "Derivative matrix: ", derivative_matrix

        # Now build the normal equations for least squares fitting
        normal_matrix = empty(shape=(numFactors, numFactors + 1))
        for i in range(0, numFactors):
            for j in range(0, numFactors):
                temp = derivative_matrix[0][i] * derivative_matrix[0][j]
                for k in range(1, numPoints):
                    temp += derivative_matrix[k][i] * derivative_matrix[k][j]

                normal_matrix[i][j] = temp

            temp = derivative_matrix[0][i] * -(zBedProbePoints[0] + corrections[0])
            for k in range(1, numPoints):
                temp += derivative_matrix[k][i] * -(zBedProbePoints[k] + corrections[k])

            normal_matrix[i][numFactors] = temp

        # print "Normal matrix: ", normal_matrix

        solution = empty(numFactors)
        gaussjordan(normal_matrix, solution, numFactors)

        for i in range(0, numFactors):
            # print solution[i]
            if isnan(solution[i]):
                raise Exception(
                    "Unable to calculate corrections. Please make sure the bed probe points are all distinct.")

        # if (debug):
        # DebugPrint(PrintVector("Solution", solution))

        # Calculate and display the residuals
        residuals = []
        for i in range(0, numPoints):
            r = zBedProbePoints[i]
            for j in range(0, numFactors):
                r += solution[j] * derivative_matrix[i][j]
                residuals.append(r)

        # DebugPrint(PrintVector("Residuals", residuals))

        deltapar.adjust(numFactors, solution, normalise)

        # Calculate the expected probe heights using the new parameters

        expectedResiduals = zeros(numPoints)
        sumOfSquares = 0.0
        for i in range(0, numPoints):
            for axis in range(0, 3):
                probe_motor_positions[i][axis] += solution[axis]

            newZ = deltapar.inverse_transform(probe_motor_positions[i][0], probe_motor_positions[i][1],
                                              probe_motor_positions[i][2])

            corrections[i] = newZ
            expectedResiduals[i] = zBedProbePoints[i] + newZ
            sumOfSquares += expectedResiduals[i] * expectedResiduals[i]

        expectedRmsError = sqrt(sumOfSquares / numPoints)
        # DebugPrint(PrintVector("Expected probe error", expectedResiduals))

        # Decide whether to do another iteration Two is slightly better than one, but three doesn't improve things.
        # Alternatively, we could stop when the expected RMS error is only slightly worse than the RMS of the residuals.
        print "Calibrated %d factors using %d points, deviation before: %f , after: %f" % (numFactors, numPoints, sqrt(
            initial_sum_of_squares / numPoints), expectedRmsError)

        iteration += 1
        if iteration == 2:
            break


def calc_probe_points(num_points=10, bed_radius=85.0):
    """Generate a set of probe points based on the printing radius.
    :param num_points: how many points
    :return: a pair of (Numpy) arrays respectively holding x and y coordinates
    """
    x = zeros(num_points)
    y = zeros(num_points)

    if num_points >= 7:
        for i in range(0, 6):
            x[i] = bed_radius * sin((2 * pi * i) / 6)
            y[i] = bed_radius * cos((2 * pi * i) / 6)

    if num_points >= 10:
        for i in range(6, 9):
            x[i] = bed_radius / 2 * sin((2 * pi * (i - 6)) / 3)
            y[i] = bed_radius / 2 * cos((2 * pi * (i - 6)) / 3)

        x[9] = 0.0
        y[9] = 0.0

    else:
        x[9] = 0.0
        y[9] = 0.0

    return x, y


def generate_commands(bot, bed_radius, firmware):
    m665 = "M665 R%.2f L%.2f" % (bot.delta_radius, bot.diagonal)
    m666 = "M666 X%.2f Y%.2f Z%.2f" % (bot.xstop, bot.ystop, bot.zstop)

    if firmware == 'RRFW':
        m665 += " H%.2f B%.2f X%.2f Y%.2f Z%.2f" % (bot.homedHeight, bed_radius, bot.xadj, bot.yadj, bot.zadj)

    elif firmware == 'Marlin':
        pass

    elif firmware == 'MarlinRC':
        m666 += " H%.2f A%.2f B%.2f C%.2f" % (bot.homedHeight, bot.xadj, bot.yadj, bot.zadj)

    elif firmware == 'Repetier':
        pass

    elif firmware == 'Smoothieware':
        m665 += " D%.2f E%.2f H%.2f Z%.2f" % (bot.xadj, bot.yadj, bot.zadj, bot.homedHeight)

    return m665, m666


def main(parser):
    args = vars(parser.parse_args())
    config_str = "Your initial setup:\n"
    for k, v in args.items():
        config_str += "%s : %s\n" % (k, v)
    print "%s\n%s" % (parser.description, config_str)

    if args['end_stops']:
        xs, ys, zs = args['end_stops'].split(',')
    else:  # default end stop setup
        xs, ys, zs = 0, 0, 0

    if args['tower_pos']:
        xt, yt, zt = args['tower_pos'].split(',')
    else:  # default tower angular positions
        xt, yt, zt = 0, 0, 0

    bot = DeltaParameters(diagonal=args['diagonal_rod_length'], delta_radius=args['delta_radius'],
                          height=args['homed_height'], xstop=float(xs), ystop=float(ys), zstop=float(zs),
                          xadj=float(xt), yadj=float(yt), zadj=float(zt))
    x, y = calc_probe_points(num_points=args['num_points'], bed_radius=args['bed_radius'])

    if args['z_probe_points']:
        z = array(args['z_probe_points'].split(','), dtype=float)
        if z.size != args['num_points']:
            print "Warning: the z probe points array does not correspond to the %d probe points: check the size!" % \
                  args['num_points']
    else:
        z = zeros(args['num_points'])

    do_delta_calibration(bot, args['num_factors'], args['num_points'], x, y, z, args['normalize'])

    m665, m666 = generate_commands(bot, args['bed_radius'], args['firmware'])
    print "Use the following G-code commands to configure your printer firmware (%s):" % args['firmware']
    print "%s\n%s\n" % (m665, m666)


if __name__ == '__main__':
    parser = ArgumentParser(
        description="""Delta 3D printer calibration utility. Based on David Crocker's (DC42) javascript code.""")

    parser.add_argument("-n", "--normalize", action="store_true", dest="normalize",
                        help="enable normalization")

    parser.add_argument('-p', "--probe_points", dest='num_points', default=10, type=int,
                        help="how many probe points are used for calibration. Default: 10")

    parser.add_argument('-f', "--cal_factors", dest='num_factors', default=6, type=int,
                        help="how many factors to calibrate. Default: 6")

    parser.add_argument('-z', "--z_probe_points", dest='z_probe_points', default=None,
                        help="A comma separated list of floats, representing the z coordinate of each probe point. \
                        If the 1st number is negative, please add a white space in front of it, such as: -z ' -3.5,8.23,...'. \
                        No whitespaces are allowed. Default: a zeros array.")

    parser.add_argument('-b', "--bed_radius", dest='bed_radius', default=85.0, type=float,
                        help="printing bed radius. Default: 85mm (Kossel mini)")

    parser.add_argument('-r', "--diagonal_rod_length", dest='diagonal_rod_length', default=215.0, type=float,
                        help="Diagonal rod length. Default: 215.0mm (Kossel mini)")

    parser.add_argument('-d', "--delta_radius", dest='delta_radius', default=105.6, type=float,
                        help="Delta radius. Default: 105.6mm (Kossel mini)")

    parser.add_argument('-H', "--homed_height", dest='homed_height', default=250.0, type=float,
                        help="Effector height when homed. Default: 105.6mm (Kossel mini)")

    parser.add_argument('-e', "--end_stops", dest='end_stops', default=None,
                        help="Height of the end stops (x,y,z) when homed. It is a comma separated list of floats, \
                             no whitespaces are allowed. Default: a zeros array")

    parser.add_argument('-t', "--tower_pos", dest='tower_pos', default=None,
                        help="Angular tower corrections (x,y,z). It is a comma separated list of floats. \
                             If the 1st number is negative, please add a white space in front of it, such as: -z ' -3.5,8.23,...'.\
                             No whitespaces are allowed. Default: a zeros array")

    parser.add_argument('-F', "--firmware", dest='firmware', default='RRFW', type=str,
                        help="Adopted firmware. Supported firmwares are: %s. Default: RRFW" % SUPPORTED_FIRMWARES)

    main(parser)
