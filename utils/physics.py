""" Physics utilities
"""

import exceptions
import math

class MeasuredQuantity(object):
    """ Class that represents a measured quantiy with errors and units
    """
    def __init__(self, value, error, units=""):
        if error < 0:
            raise exceptions.ValueError()
        if math.isnan(value) or math.isnan(error):
            raise exceptions.ValueError()
        if units == "%":
            value *= 100
            error *= 100
        self.value = value
        self.error = error
        self.units = units
        self.force_exponent = None

    def is_significant(self):
        """ Checks if the uncertainty is smaller than the value
        """
        if self.value == 0:
            return self.error == 0
        return self.error / math.fabs(self.value) < 1.0

    def calculate_precision(self):
        """ Calculates the number of significant digits from the uncertainty
        """
        if self.error == 0:
            return float("inf")

        err_log_10 = math.log10(self.error)
        if err_log_10 < 0:
            return int(math.floor(err_log_10))
        else:
            return int(math.ceil(err_log_10))

    def to_string(self, plus_minus="#pm"):
        """ String representation
        """
        precision = self.calculate_precision()
        if math.isinf(precision):
            precision = 6
        abs_precision = abs(precision)

        if self.force_exponent is not None:
            exp = self.force_exponent
        elif abs_precision > 3:
            if precision > 0:
                exp = precision - 1
            else:
                exp = precision + 1
        else:
            exp = None

        if exp is not None:
            value = self.value / (10**exp)
            error = self.error / (10**exp)
            format_string = "({{value:.{prec}f}} {{pm}} {{error:.{prec}f}}) #times 10^{{{{{{exp}}}}}}".\
                format(prec=abs(precision-exp))
            result = format_string.\
                format(value=value, error=error, exp=exp, pm=plus_minus)
        else:
            format_string = "{{value:.{prec}f}} {{pm}} {{error:.{prec}f}}".\
                format(prec=abs_precision)
            result = format_string.\
                format(value=self.value, error=self.error,pm=plus_minus)
        if self.units:
            result += " {}".format(self.units)
        return result

    def __str__(self):
        return self.to_string("+/-")

    def __add__(self, other):
        if isinstance(other, MeasuredQuantity):
            if self.units != other.units:
                return NotImplemented
            return MeasuredQuantity(self.value + other.value, math.sqrt(self.error**2 + other.error**2), self.units)
        elif isinstance(other, (int, float)):
            if self.units:
                return NotImplemented
            return MeasuredQuantity(self.value + other, self.error)
        return NotImplemented

    def __neg__(self):
        return MeasuredQuantity(-self.value, self.error, self.units)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if isinstance(other, MeasuredQuantity):
            value = self.value * other.value
            error = math.sqrt(self.error**2 / self.value**2 + other.error**2 / other.value**2) * value
            units = self.units
            if units and other.units:
                if units == other.units:
                    if units == "%":
                        value /= 100
                        error /= 100
                    else:
                        units = "({})^2".format(units)
                else:
                    units += " " + other.units
            return MeasuredQuantity(value, error, units)
        elif isinstance(other, (int, float)):
            return MeasuredQuantity(self.value * other, self.error * other)
        else:
            return NotImplemented

    def __invert__(self):
        if self.units != "%":
            return MeasuredQuantity(1 / self.value, self.error / self.value**2, self.units)
        return NotImplemented

    def __div__(self, other):
        return self * (~other)

    def __lt__(self, other):
        if isinstance(other, MeasuredQuantity):
            return self.value + self.error < other.value - other.error
        elif isinstance(other, (int, float)):
            return self.value + self.error < other
        return NotImplemented

    def __gt__(self, other):
        if isinstance(other, MeasuredQuantity):
            return self.value + self.error > other.value - other.error
        elif isinstance(other, (int, float)):
            return self.value + self.error > other
        return NotImplemented

    def ___le__(self, other):
        return not self > other

    def __ge__(self, other):
        return not self < other

    def __ne__(self, other):
        return self > other or self < other

    def __eq__(self, other):
        return not self != other
