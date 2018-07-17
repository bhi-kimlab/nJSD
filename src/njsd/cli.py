"""
Module that contains the command line app.

Why does this file exist, and why not put this in __main__?

  You might be tempted to import things from __main__ later, but that will cause
  problems: the code will get executed twice:

  - When you run `python -mnjsd` python will execute
    ``__main__.py`` as a script. That means there won't be any
    ``njsd.__main__`` in ``sys.modules``.
  - When you import __main__ it will get executed again (as a module) because
    there's no ``njsd.__main__`` in ``sys.modules``.

  Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""
import argparse
import njsd

# Create top-level parser.
parser = argparse.ArgumentParser(description='Calculate network-based Jensen-Shannon Divergence.')
parser.add_argument('-n', '--network', required=True, help='Pre-defined network')
parser.add_argument('-r', '--ref', required=True, help='Reference gene expression profile')
parser.add_argument('-q', '--query', required=True, help='Query gene expression profile')
parser.add_argument('-o', '--output', required=True, help='Output file.')
parser.add_argument('-t', '--geneset', default=None, help='Gene set list')


def main(args=None):
    args = parser.parse_args(args=args)
    if args.geneset is not None:
        # Run gene set-specified nJSD calculation.
        result = njsd.njsd_geneset(args.network, args.ref, args.query, args.geneset, file=args.output)
    else:
        # Run transcriptome-wide nJSD calculation.
        result = njsd.njsd_all(args.network, args.ref, args.query, file=args.output)
