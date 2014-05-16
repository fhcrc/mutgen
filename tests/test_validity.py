from helpers import TestGenerator

import random
random.seed("I can get you a toe dude. With nail polish!")


def mut_pattern(motif, prob):
    # XXX - need to replace this with something more intelligent
    return [motif, 0, [('A', prob), ('C', 0), ('G', 0), ('T', 0)]]


test_generator = TestGenerator(sample_size = 1000, almost_equal_places = 1)

TestSimpleSinglePattern = test_generator.make_test("TestSimpleSinglePattern",
    'CCCGTCCCCGACCCCCGGCCCCC',
    [mut_pattern('GT', 0.2),
        mut_pattern('GA', 0.3),
        mut_pattern('GG', 0.7)],
    dict(GT=0.2, GA=0.3, GG=0.7))

