import os
import sys

import pypsa

if __name__ == "__main__":
    # Read file name as first command line argument
    fn_in = sys.argv[1]
    fn_out = sys.argv[2]
    # If fn_out already exists, do nothing
    if not os.path.exists(fn_out):
        n = pypsa.Network(fn_in)
        # Make sure the directory containing the output file exists
        os.makedirs(os.path.dirname(fn_out), exist_ok=True)
        print(f"Compressing {fn_in} to {fn_out}")
        n.export_to_netcdf(
            fn_out, compression=dict(zlib=True, complevel=3, least_significant_digit=3)
        )
