from helpers import TestGenerator

import random
random.seed("I can get you a toe dude. With nail polish!")

test_generator = TestGenerator(sample_size = 1000, almost_equal_places = 1)


TestSimpleSinglePattern = test_generator.make_test("TestSimpleSinglePattern",
    # For the frequency math to work, each kmer result tested for should only show up once here
    # (this requirement could be fixed...)
    'CCCGTCCCCGACCCCCGGCCCCC',
    # Mutation model specification
    [   ['GT', 0, dict(A=0.2)],
        ['GA', 0, dict(A=0.3)],
        ['GG', 0, dict(A=0.7)]],
    # Specification of expected empirical frequencies
    dict(GT=0.2, GA=0.3, GG=0.7))


TestDifferentMutTargets = test_generator.make_test("TestDifferentMutTargets",
    'CCCGTCCCCGACCCCCGGCCCCC',
    [   ['GT', 0, dict(A=0.2, C=0.5)],
        ['GA', 0, dict(A=0.3, C=0.2)],
        ['GG', 0, dict(A=0.7, C=0.1)]],
    # Specification of expected empirical frequencies; this time with mutated_to nt specification
    {   ('GT', 'A'): 0.2, ('GT', 'C'): 0.5,
        ('GA', 'A'): 0.3, ('GA', 'C'): 0.2,
        ('GG', 'A'): 0.7, ('GG', 'C'): 0.1})


TestSimpleAmbiguous = test_generator.make_test("TestSimpleAmbiguous",
    'CCCGTCCCCGACCCCCGGCCCCC',
    [   ['GW', 0, dict(A=0.2)],
        ['GG', 0, dict(A=0.7)]],
    dict(GT=0.2, GA=0.2, GG=0.7))


# We add padding to each of the motif patterns so that the mut_index positions line up. Want to make sure this
# is working properly; Also tests ambiguous a bit.
TestJaggedMotifs = test_generator.make_test("TestSimpleAmbiguous",
    'CCCGTCCCCGACCCCCGGCCCCC',
    [   ['CGW', 1, dict(A=0.2)],
        ['GGS', 0, dict(A=0.7)]],
    dict(CGTC=0.2, CGAC=0.2, CGGC=0.7))


