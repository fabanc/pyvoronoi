import time, logging
import pyvoronoi

def populate_segments(max_x, max_y):
    pv = pyvoronoi.Pyvoronoi()
    count_segment = 0
    for x in xrange(max_x):
        for y in xrange(max_y):
            pv.AddSegment([[x,y], [x, y + 1]])
            count_segment += 1

    time_before = time.time()
    pv.Construct()
    time_after = time.time()
    logging.info("Run pyvoronoi. Time (sec): {0}. Number of input points: {1} - segments: {2}".format(
            time_after - time_before,
            0,
            count_segment
    ))

    logging.info("Count output structures. Vertices: {0}, Edges: {1}, Cells: {2}".format(
        pv.CountVertices(),
        pv.CountEdges(),
        pv.CountCells(),
    ))

    del pv

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s - Script: %(filename)s - Line: %(lineno)d - %(levelname)s - %(message)s',
                        datefmt='%m-%d %H:%M:%S')
    populate_segments(100, 12000)
    #populate_segments(100, 10000)
