from helpers import TestGenerator

import random
random.seed("I can get you a toe dude. With nail polish!")

test_generator = TestGenerator(sample_size = 1000, almost_equal_places = 1)


TestSimpleSinglePattern = test_generator.make_test("TestSimpleSinglePattern",
    'CCCGTCCCCGACCCCCGGCCCCC',
    # Mutation model specification
    [   ['GT', 0, dict(A=0.2)],
        ['GA', 0, dict(A=0.3)],
        ['GG', 0, dict(A=0.7)]],
    # Specification of expected empirical frequencies
    dict(GT=0.2, GA=0.3, GG=0.7))


TestDifferentMutTargets = test_generator.make_test("TestDifferentMutTargets",
    'CCCGTCCCCGACCCCCGGCCCCC',
    # Mutation model specification
    [   ['GT', 0, dict(A=0.2, C=0.5)],
        ['GA', 0, dict(A=0.3, C=0.2)],
        ['GG', 0, dict(A=0.7, C=0.1)]],
    # Specification of expected empirical frequencies
    {   ('GT', 'A'): 0.2, ('GT', 'C'): 0.5,
        ('GA', 'A'): 0.3, ('GA', 'C'): 0.2,
        ('GG', 'A'): 0.7, ('GG', 'C'): 0.1})

