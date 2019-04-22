proc color_mononucleosome {} {
    mol delrep 0 top

    # DNA
    mol color ColorID 8
    mol representation NewCartoon 0.300000 10.000000 4.100000 0
    mol selection (residue 974 to 1268)
    mol material Opaque
    mol addrep top

    # H3
    mol color ColorID 0
    mol representation NewCartoon 0.300000 10.000000 4.100000 0
    mol selection (residue 0 to 134 or residue 487 to 621)
    mol material Opaque
    mol addrep top

    # H4
    mol color ColorID 1
    mol representation NewCartoon 0.300000 10.000000 4.100000 0
    mol selection (residue 135 to 236 or residue 622 to 723)
    mol material Opaque
    mol addrep top

    # H2A
    mol color ColorID 4
    mol representation NewCartoon 0.300000 10.000000 4.100000 0
    mol selection (residue 237 to 364 or residue 724 to 851)
    mol material Opaque
    mol addrep top

    # H2B
    mol color ColorID 3
    mol representation NewCartoon 0.300000 10.000000 4.100000 0
    mol selection (residue 365 to 486 or residue 852 to 973)
    mol material Opaque
    mol addrep top
}