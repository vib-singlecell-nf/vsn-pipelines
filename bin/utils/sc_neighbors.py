def add_options(parser, method):
    parser.add_option(
        "-n", "--n-neighbors",
        type="int",
        action="store",
        dest="n_neighbors",
        default=15,
        help="[{}], The size of local neighborhood (in terms of number of neighboring data points) used for manifold"
             " approximation.".format(method)
    )
    parser.add_option(
        "-p", "--n-pcs",
        type="int",
        action="store",
        dest="n_pcs",
        default=30,
        help="[{}], Use this many PCs.".format(method)
    )
    return parser
