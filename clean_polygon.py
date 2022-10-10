def clean_polygon(fn):
    """ Clean WKT """
    from shapely.wkt import loads, dumps
    p = loads(open(fn, "r").read())
    p = p.buffer(0)
    with open(fn, "w") as tf:
        tf.write(dumps(p))


