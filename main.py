from math import ceil, log2


def remove_trailing_zeroes(ls):
    while (len(ls) > 1 and ls[-1] == 0):
        ls.pop()
    return ls


def FFT(polynomial, number_of_values, modulo, base_w, base_n):
    def evaluate():
        auxiliary_w = [None] * number_of_values
        val = 1
        for i in range(number_of_values):
            auxiliary_w[i] = val
            val = (w * val) % modulo
        polynomial_evaluation = polynomial[::]
        count = 2
        while (count <= number_of_values):
            dist = number_of_values // count
            new = polynomial_evaluation[::]
            for start in range(dist):
                for i in range(count // 2):
                    position = 2 * i * dist + start
                    a = polynomial_evaluation[position]
                    b = polynomial_evaluation[position + dist]
                    new[i * dist + start] = (a + auxiliary_w[i * dist] * b) % modulo
                    new[(i + count // 2) * dist + start] = (a - auxiliary_w[i * dist] * b) % modulo
            count *= 2
            polynomial_evaluation = new
        return polynomial_evaluation

    number_of_values = 2 ** ceil(log2(number_of_values))
    polynomial.extend([0] * (number_of_values - len(polynomial)))
    w = pow(base_w, base_n // number_of_values, modulo)
    return evaluate()


def multiply_polynomials(a, b, modulo=3 * 2 ** 30 + 1, base_w=125, base_n=2 ** 30):
    l = len(a) + len(b) - 1
    evaluated_a = FFT(a, l, modulo, base_w, base_n)
    evaluated_b = FFT(b, l, modulo, base_w, base_n)
    evaluated_product = []
    for i in range(len(evaluated_a)):
        evaluated_product.append((evaluated_a[i] * evaluated_b[i]) % modulo)
    product = FFT(evaluated_product, len(evaluated_product), modulo, pow(base_w, -1, modulo), base_n)
    d = pow(len(evaluated_product), -1, modulo)
    for i in range(len(product)):
        product[i] = (product[i] * d) % modulo
    return remove_trailing_zeroes(product)
