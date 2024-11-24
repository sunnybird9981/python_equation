#########################################################
# Class representing a monomial
#########################################################
class Monomial():
    def __init__(self, coefficient, variables : dict):
        self.coefficient = coefficient
        self.variables = variables
        self.vars_list = self.get_variables_list()
        if self.coefficient == 0:
            self.variables = {}

    # assign value for var
    # values arg must be {'var1': 'value1', 'var2': 'value2', ...}
    def evaluate(self, values):
        result = self.copy()

        for var, value in values.items():
            if var in result.variables:
                result *= (value ** self.variables.get(var))
                if result.coefficient == 0:
                    break
                else:
                    del result.variables[var]

        return result

    # partial_derivative a monomial
    def partial_derivative(self, var):
        if (var not in self.variables) or self.variables[var] == 0:
            return Monomial(0, {})
        else:
            new_coefficient = self.coefficient * self.variables[var]
            new_variables = self.variables.copy()
            if new_variables[var] == 1:
                del new_variables[var]
            else:
                new_variables[var] -= 1

            return Monomial(new_coefficient, new_variables)

    # return chars set used as variables in this monomial
    def get_variables_set(self):
        variables = set()

        variables.update(self.variables.keys())

        return variables

    # return sorted chars list used as variables in this monomial
    def get_variables_list(self):
        return sorted(self.get_variables_set())

    # check if self and "other" monomials are like terms
    def is_like_terms(self, other):
        if self.variables == other.variables:
            return True
        else:
            return False

    def delete_variables_if_constant(self):
        if self.is_constant():
            self.variables.clear()

    # check if self contains no variables or not
    def is_constant(self):
        for key, value in self.variables.items():
            if value != 0:
                return False
        else:
            return True

    # if this monomial is 0, set this monomial to Monomial(0, {})
    def delete_variables_if_zero(self):
        if self.is_zero():
            self.variables.clear()

    # check if coefficient is 0 or not
    def is_zero(self):
        if self.coefficient == 0:
            return True
        else:
            return False

    # shallow copy
    def copy(self):
        return Monomial(self.coefficient, self.variables.copy())

    # evaluate exponentiation
    def mul_power(self, mult):
        new_monomial = self.copy()

        for var, power in new_monomial.variables.items():
            new_monomial.variables[var] *= mult

        return new_monomial


    # add a monomial or a polynomial
    def add(self, other):
        if isinstance(other, Monomial):
            return self.add_monomial(other)
        elif isinstance(other, Polynomial):
            return other.add_monomial(self)
        elif isinstance(other, (int, float, complex)):
            return self.add_constant(other)
        else:
            raise ValueError('add(other) : "other" must be Polynomail or Monomial.')

    # add a monomial
    def add_monomial(self, other):
        if self.is_like_terms(other):
            return self.add_like_terms(other)
        else:
            return self.add_unlike_terms(other)

    # add a monomial which and self are like terms
    def add_like_terms(self, other):
        new_monomial = Monomial(self.coefficient, self.variables.copy())
        new_monomial.coefficient += other.coefficient
        new_monomial.delete_variables_if_zero()

        return Polynomial([new_monomial])

    # add a monomial which and self are not like terms
    def add_unlike_terms(self, other):
        new_monomials = [Monomial(self.coefficient, self.variables.copy())]
        new_monomials.append(other)

        return Polynomial(new_monomials)

    # add a polynomial
    def add_polynomial(self, other):
        if isinstance(other, Polynomial):
            return other.add_monomial(self)
        else:
            raise ValueError('add_polynomial(other) : "other" must be Polynomail.')

    # add a number which contains no variables
    def add_constant(self, other):
        other_as_monomial = Monomial(other, {})

        return self.add_monomial(other_as_monomial)

    # subtract a monomial or a polynomial
    def sub(self, other):
        if isinstance(other, Monomial):
            return self.add_monomial(- other)
        elif isinstance(other, Polynomial):
            return self.add_polynomial(-other)
        elif isinstance(other, (int, float, complex)):
            return self.add_constant(-other)
        else:
            raise ValueError('sub(other) : "other" must be Polynomail or Monomial.')

    # multiply a number
    def mul(self, other):
        if isinstance(other, (int, float, complex)):
            return self.mul_constant(other)
        elif isinstance(other, Monomial):
            return self.mul_monomial(other)
        elif isinstance(other, Polynomial):
            return self.mul_polynomial(other)
        else:
            raise ValueError('mul(other) : "other" must be int, float, Monomial or Polynomial')

    # multiply by a constant
    def mul_constant(self, other):
        new_coefficient = self.coefficient * other

        return Monomial(new_coefficient, self.variables)

    # multiply a Monomial
    def mul_monomial(self, other):
        new_coefficient = self.coefficient * other.coefficient
        new_variables = self.variables.copy()

        for var, power in other.variables.items():
            if var in new_variables:
                if new_variables[var] == - power:
                    del new_variables[var]
                else:
                    new_variables[var] += power
            else:
                new_variables[var] = power

        return Monomial(new_coefficient, new_variables)

    # multiply a Polynomial
    def mul_polynomial(self, other):
        new_polynomial = Polynomial([])

        for m in other.monomials:
            new_polynomial.monomials.append(self.mul_monomial(m))

        new_polynomial.combine_like_terms()

        return new_polynomial

    # evaluate -(self)
    def rsub(self, other):
        minus_self = -self

        return minus_self.add(other)

    # evaluate exponentiation
    def pow_constant(self, other):
        new_monomial = self.mul_power(other)
        new_monomial.coefficient **= other

        return new_monomial

    # evaluate exponentiation
    def pow(self, other):
        if isinstance(other, (int, float)):
            return self.pow_constant(other)
        else:
            raise ValueError('power must be int or float')

    def __add__(self, other):
        return self.add(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.sub(other)

    def __rsub__(self, other):
        return self.rsub(other)

    def __mul__(self, other):
        return self.mul(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, other):
        return self.pow(other)

    def __neg__(self):
        negated_monomial = Monomial(- self.coefficient, self.variables.copy())

        return negated_monomial

    def __hash__(self):
        return hash((self.coefficient, frozenset(self.variables.items())))

    # decide the format for print()
    # for example, print(Monomial(3, {"x" : 2})) return "3x^2"
    def __repr__(self):
        vars = self.variables.items()
        formatted_str = ""

        if any(vars):
            if self.coefficient == -1:
                formatted_str += "-"
            elif self.coefficient != 1:
                formatted_str += str(self.coefficient)

            for var, power in self.variables.items():
                if power != 1:
                    formatted_str += f"{var}^{power}"
                else:
                    formatted_str += var
        else:
            formatted_str = str(self.coefficient)

        return formatted_str

#########################################################
# Class representing a polynomial
#########################################################
class Polynomial():
    def __init__(self, monomials : list):
        self.monomials = monomials
        self.vars_list = self.get_variables_list()

    # assign value for var
    # values arg must be {'var1': 'value1', 'var2': 'value2', ...}
    def evaluate(self, values):
        new_polynomial = Polynomial([])

        for m in self.monomials:
            new_polynomial.monomials.append(m.evaluate(values))

        new_polynomial.combine_like_terms()

        return new_polynomial

    # partial_derivative this polynomial
    def partial_derivative(self, var):
        new_polynomial = []

        for m in self.monomials:
            m_partial_derivative = m.partial_derivative(var)

            if m_partial_derivative.coefficient != 0:
                new_polynomial.append(m_partial_derivative)

        if new_polynomial:
            return Polynomial(new_polynomial)
        else :
            return Polynomial([Monomial(0, {})])

    # shallow copy
    def copy(self):
        new_monomials = self.monomials[:]

        return Polynomial(new_monomials)

    # count max degree of var in this polynomial
    def get_degree_in_variable(self, var):
        max_degree = 0

        for m in self.monomials:
            max_degree = max(max_degree, m.variables.get(var, 0))

        return max_degree

    # sort each monomials by ascending order with respect to "var"
    def sort_by_ascending_order(self, var):
        self.monomials.sort(key=lambda m: m.variables.get(var, 0))

    # sort each monomials by descending order with respect to "var"
    def sort_by_descending_order(self, var):
        self.monomials.sort(key=lambda m: m.variables.get(var, 0), reverse=not False)

    # get chars set used as variations in this polynomial
    def get_variables_set(self):
        variables = set()

        for m in self.monomials:
            variables.update(set(m.vars_list))

        return variables

    # return sorted chars list used as variables in this monomial
    def get_variables_list(self):
        return sorted(self.get_variables_set())

    # add a polynomial or a monomial
    def add(self, other):
        if isinstance(other, Monomial):
            return self.add_monomial(other)
        elif isinstance(other, Polynomial):
            return self.add_polynomial(other)
        elif isinstance(other, (int, float, complex)):
            return self.add_constant(other)
        else:
            raise ValueError('only polynomials or monomials objects can be added to a polynomials object.')

    # add a constant, not a variable
    def add_constant(self, other):
        other_as_monomial = Monomial(other, {})

        return self.add_monomial(other_as_monomial)

    # add a monomial
    def add_monomial(self, other: Monomial):
        new_monomials = self.monomials[:]

        for i, m in enumerate(self.monomials):
            if m.is_like_terms(other):
                new_monomials[i] = m.add_like_terms(other).monomials[0]
                break
        else:
            new_monomials.append(other)

        return Polynomial(new_monomials)

    # add a monomial which and all of monomials in self are not like terms
    def add_non_matching_monomial(self, other: Monomial):
        new_monomials = self.monomials[:]
        new_monomials.append(other)

        return Polynomial(new_monomials)

    # add a polynomial
    def add_polynomial(self, other):
        new_monomials = self.monomials[:]

        for n in other.monomials:
            for i, m in enumerate(self.monomials[:]):
                if m.is_like_terms(n):
                    new_monomials[i] = m.add_like_terms(n).monomials[0]
                    break
            else:
                new_monomials.append(n)

        return Polynomial(new_monomials)

    # subtract a polynomial
    def sub_polynomial(self, other):
        new_monomials = self.monomials[:]

        for n in other.monomials:
            for i, m in enumerate(self.monomials[:]):
                if m.is_like_terms(n):
                    new_monomials[i] = m.add_like_terms(-n).monomials[0]
                    break
            else:
                new_monomials.append(-n)

        return Polynomial(new_monomials)

    # subtract a monomial or polynomial
    def sub(self, other):
        if isinstance(other, Monomial):
            return self.add_monomial(-other)
        elif isinstance(other, Polynomial):
            return self.sub_polynomial(other)
        elif isinstance(other, (int, float, complex)):
            return self.add_constant(-other)
        else:
            raise ValueError('only polynomials or monomials objects can be subtracted to a polynomials object.')


# if this polynomial contains monomials which are like terms,
# add them together.
    def combine_like_terms(self):
        term_dict = {}
        new_monomials = []

        for m in self.monomials:
            var_tuple = tuple(m.variables.items())

            if var_tuple in term_dict:
                term_dict[var_tuple] += m.coefficient
            else:
                term_dict[var_tuple] = m.coefficient

        for variables, coefficient in term_dict.items():
            if coefficient != 0:
                var_dict = dict(variables)
                new_monomials.append(Monomial(coefficient, var_dict))

        self.monomials = new_monomials

    # get term whose variables are the same as vars in "dict"
    def get_term_by_vars(self, var : dict):
        for m in self.monomials:
            if m.variables == var:
                break
        else:
            return Monomial(0, {})

        return m

    # get coefficients list
    # x^4 + 2x^2 - x + 1
    # [4, 0, 2, -1, 1]
    def get_coefficients_list(self, var):
        degree = self.get_degree_in_variable(var)
        coefficients = [0 for _ in range(degree + 1)]

        for m in self.monomials:
            index = degree - m.variables.get(var, 0)
            coefficients[index] = m.coefficient

        return coefficients

    # multiply by a number
    def mul(self, other):
        if isinstance(other, (int, float, complex)):
            return self.mul_constant(other)
        elif isinstance(other, Monomial):
            return self.mul_monomial(other)
        elif isinstance(other, Polynomial):
            return self.mul_polynomial(other)
        else:
            raise ValueError('mul(other) : "other" must be int, float, Monomial or Polynomial')

    # multiply by a constant, not a variable
    def mul_constant(self, other):
        new_polynomial = []

        for m in self.monomials:
            new_polynomial.append(m.mul_constant(other))

        return Polynomial(new_polynomial)

    # multiply by a Monomial
    def mul_monomial(self, other):
        new_polynomial = []

        for m in self.monomials:
            new_polynomial.append(m.mul_monomial(other))

        return Polynomial(new_polynomial)

    # multiply by a Polynomial
    def mul_polynomial(self, other):
        new_polynomial = Polynomial([])

        for m in self.monomials:
            m_mul_other = m.mul_polynomial(other)
            new_polynomial = new_polynomial.add_polynomial(m_mul_other)

        return new_polynomial

    # evaluate (self)^(other) : other is int
    def pow_int(self, other):
        if other > 0:
            return self.pow_natural_number(other)
        elif other == 0:
            new_polynomial = Polynomial([
                Monomial(0, {})
                ])

            return new_polynomial
        else:
            raise ValueError('other must be equal or greater than 0')

    # evaluate (self)^(other) : other is natural number
    def pow_natural_number(self, other):
        new_polynomial = self.copy()

        for i in range(other - 1):
            new_polynomial = new_polynomial.mul_polynomial(self)

        return new_polynomial

    # evaluate (self)^(other)
    def pow(self, other):
        if isinstance(other, int):
            return self.pow_int(other)
        else:
            raise ValueError('power must be int')

    # evaluate -(self)
    def rsub(self, other):
        minus_self = -self

        return minus_self.add(other)

    def __add__(self, other):
        return self.add(other)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        return self.sub(other)

    def __rsub__(self, other):
        return self.rsub(other)

    def __mul__(self, other):
        return self.mul(other)

    def __pow__(self, other):
        return self.pow(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __neg__(self):
        negated_monomials = []

        for m in self.monomials:
            dict_copy = {key: value for key, value in m.variables.items()}
            negated_monomials.append(Monomial(- m.coefficient, dict_copy))

        return Polynomial(negated_monomials)

    # decide the format for print()
    def __repr__(self):
        formatted_str = []
        for m in self.monomials:
            if formatted_str:
                if type(m.coefficient) is complex:
                    formatted_str.append(" + " + str(m))
                else:
                    if m.coefficient < 0:
                        formatted_str.append(" - " + str(m)[1:])
                    else:
                        formatted_str.append(" + " + str(m))
            else:
                formatted_str.append(str(m))
        return "".join(formatted_str)

#########################################################
# Class representing a equation
#########################################################
class Equation():
    def __init__(self, pl, pr):
        self.lhs = self._cast_to_polynomial(pl)
        self.rhs = self._cast_to_polynomial(pr)
        self.vars_list = self.get_variables_list()

    # cast self.lhs and self.rhs to Polynomial Class
    def _cast_to_polynomial(self, value):
        if isinstance(value, Polynomial):
            return value
        elif isinstance(value, Monomial):
            return Polynomial([value])
        elif isinstance(value, (int, float, complex)):
            return Polynomial([Monomial(value, {})])
        else:
            raise TypeError('unsupported type for equation side: {}'.format(type(value)))

    # shallow copy
    def copy(self):
        return Equation(self.lhs.copy(), self.rhs.copy())

    # add the same value to each hand side
    def add(self, value):
        new_lhs = self.lhs + value
        new_rhs = self.rhs + value

        return Equation(new_lhs, new_rhs)

    # subtract value from both sides
    def sub(self, value):
        new_lhs = self.lhs - value
        new_rhs = self.rhs - value

        return Equation(new_lhs, new_rhs)

    # multiply both sides
    def mul(self, value):
        new_lhs = self.lhs * value
        new_rhs = self.rhs * value

        return Equation(new_lhs, new_rhs)

    # divide both sides
    def div(self, value):
        if value != 0:
            new_lhs = self.lhs * (1 / value)
            new_rhs = self.rhs * (1 / value)

            return Equation(new_lhs, new_rhs)
        else:
            raise ValueError('cannot be divided by zero.')

    # move right_hand terms to left_hand side to form eq = 0.
    def simplify(self):
        self.lhs -= self.rhs
        self.rhs = 0
        self.lhs.combine_like_terms()

    # get chars set used as variations in this equation
    def get_variables_set(self):
        variables = set()

        variables.update(set(self.lhs.vars_list))
        variables.update(set(self.rhs.vars_list))

        return variables

    # return sorted chars list used as variables in this monomial
    def get_variables_list(self):
        return sorted(self.get_variables_set())

    # get max degree of a variable in this equation
    def get_degree_in_variable(self, var):
        return max(self.lhs.get_degree_in_variable(var),
                   self.rhs.get_degree_in_variable(var))

    def solve_linear_eq_for_one_var(self, var):
        self.simplify()
        self.lhs.sort_by_descending_order(var)
        self.rhs = self.lhs.monomials[0] - self.lhs
        self.lhs = Polynomial([self.lhs.monomials[0]])
        self.rhs *= (1 / self.lhs.monomials[0].coefficient)
        self.lhs *= (1 / self.lhs.monomials[0].coefficient)
        self.rhs.combine_like_terms()

    def __repr__(self):
        return f"{self.lhs} = {self.rhs}"
