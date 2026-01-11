import math
import random
import matplotlib.pyplot as plt

def signed_area(a, b, c):
    """
    Calculates the signed area of a triangle formed by three points.
    
    Args:
        a : tuple[float, float] : First vertex of the triangle.
        b : tuple[float, float] : Second vertex of the triangle.
        c : tuple[float, float] : Third vertex of the triangle.
    
    Returns:
        float : Signed area of the triangle.
                Positive if vertices are in counterclockwise order.
                Negative if vertices are in clockwise order.
                Zero if points are collinear.
    """
    return ((b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]))/2 

def dist(p, q):
    """
    Calculates the Euclidean distance between two points.
    
    Args:
        p : tuple[float, float] : First point.
        q : tuple[float, float] : Second point.
    
    Returns:
        float : Euclidean distance between p and q.
    """
    return math.sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2)

def test_segment_intersection(p, q):
    """
    Tests if two line segments intersect.
    
    Args:
        p : tuple[tuple[float, float], tuple[float, float]] : First segment defined by two points.
        q : tuple[tuple[float, float], tuple[float, float]] : Second segment defined by two points.
    
    Returns:
        bool : True if the segments intersect (including endpoint touch).
               False otherwise.
    """
    if (signed_area(p[0], p[1], q[0])*signed_area(p[0], p[1], q[1]) < 0) and (signed_area(q[0], q[1], p[0])*signed_area(q[0], q[1], p[1]) < 0):
        return True
    else:
        if on_segment_II(p[0], q) or on_segment_II(p[1], q) or on_segment_II(q[0], p) or on_segment_II(q[1], p):
            return True
        else:
            return False
        
def on_segment_II(p, s):
    """
    Checks if a point p lies on a segment s using coordinate bounds.
    
    Args:
        p : tuple[float, float] : Point to check.
        s : tuple[tuple[float, float], tuple[float, float]] : Segment defined by two endpoints.
    
    Returns:
        bool : True if p lies on segment s.
               False otherwise.
    """
    if signed_area(p, s[0], s[1]) != 0:
        return False
    if s[0][0] <= p[0] and p[0] <= s[1][0] or s[0][0] >= p[0] and p[0] >= s[1][0]:
        if s[0][1] <= p[1] and p[1] <= s[1][1] or s[0][1] >= p[1] and p[1] >= s[1][1]:
            return True
    return False

def on_segment(p, s):
    """
    Checks if a point p lies on a segment s using distance equality.
    
    Args:
        p : tuple[float, float] : Point to check.
        s : tuple[tuple[float, float], tuple[float, float]] : Segment defined by two endpoints.
    
    Returns:
        bool : True if the sum of distances from p to both endpoints equals
               the segment length (i.e., p is on the segment).
               False otherwise.
    """
    if signed_area(p, s[0], s[1]) != 0:
        return False
    return float(dist(p, s[0])) + float(dist(p, s[1])) == float(dist(s[0], s[1]))

def in_triangle(p, t):
    """
    Determines if a point p is inside a triangle t using signed areas.
    
    Args:
        p : tuple[float, float] : Point to check.
        t : list[tuple[float, float]] : Triangle defined by three vertices.
    
    Returns:
        bool : True if p is inside the triangle (or on its boundary).
               False otherwise.
    """
    a1 = signed_area(t[0], t[1], p)
    a2 = signed_area(t[1], t[2], p)
    a3 = signed_area(t[2], t[0], p)

    if a1 >= 0 and a2 >= 0 and a3 >= 0 or a1 <= 0 and a2 <= 0 and a3 <= 0:
        return True
    else:
        return False

def y_min(p):
    """
    Finds the point with the minimum y-coordinate (and minimum x if tied).
    
    Args:
        p : list[tuple[float, float]] : List of points.
    
    Returns:
        tuple[float, float] : Point with the smallest y-coordinate.
                              In case of tie, returns the one with smallest x-coordinate.
    """
    return min(p, key=lambda x: [x[1], x[0]])

# SECTION 1

def sign(x, eps=1e-12):
    """
    Returns the sign of a real number with tolerance.

    Args:
        x : float : Number whose sign is to be determined.
        eps : float : Margin to consider x ≈ 0.

    Returns:
        int : 1 if x > eps.
              -1 if x < -eps.
               0 if |x| <= eps.
    """
    if x > eps:
        return 1
    if x < -eps:
        return -1
    return 0


def polygon_orientation(p):
    """
    Determines the orientation of a polygon.

    Args:
        p : list[tuple[float, float]] : List of polygon vertices.

    Returns:
        int : +1 if the polygon is counterclockwise.
              -1 if it is clockwise.
               0 if the area is approximately zero (degenerate).
    """
    return sign(polygon_area(p))


def is_convex(p):
    """
    Checks if a polygon is convex.

    Args:
        p : list[tuple[float, float]] : List of polygon vertices in order.

    Returns:
        bool : True if all consecutive turns have the same sign,
               ignoring null turns (collinear points).
    """
    n = len(p)
    turns = []
    for i in range(n):
        a, b, c = p[i-1], p[i], p[(i+1)%n]
        turns.append(sign(signed_area(a, b, c)))

    non_zero_turns = [g for g in turns if g != 0]
    return all(g == non_zero_turns[0] for g in non_zero_turns) if non_zero_turns else True

    
    
# SECTION 2
def section_2_point_in_polygon_method1_triangles(P, q):
    """
    Determines if a point 'q' is inside a convex polygon by
    decomposing the polygon into triangles.

    Args:
        P : list[list[float, float]]
            List of polygon vertices in order (assumed convex).
        q : list[float, float]
            Point to check.

    Returns:
        bool :
            True if 'q' is inside any of the triangles
            formed by (P[0], P[i], P[i+1]).
            False otherwise.
    """
    p0 = P[0]
    for i in range(1, len(P) - 1):
        vi_1 = P[i]
        vi_2 = P[i+1]
        triangle = [p0, vi_1, vi_2]

        if in_triangle(q, triangle):
            return True

    return False

def section_2_point_in_polygon_method2_ray_casting(polygon_points, point, p_far):
    """
    Determines if a point is inside a polygon using the
    ray casting algorithm (counting intersections with a ray).

    Args:
        polygon_points : list[list[float, float]]
            List of polygon vertices.
        point : list[float, float]
            Point to check.
        p_far : list[float, float]
            Far point that defines the ray direction.

    Returns:
        tuple :
            is_inside : bool
                True if the point is inside the polygon or on a side.
            n_cuts : int
                Number of intersections of the ray with polygon sides.
    """
    is_inside = False
    n_cuts = 0
    line = [point, p_far]

    for i in range(len(polygon_points) - 1):
        segment = [polygon_points[i], polygon_points[i+1]]
        intersection = test_segment_intersection(line, segment)

        if intersection:
            n_cuts += 1
        if on_segment(point, segment):
            is_inside = True

    segment = [polygon_points[-1], polygon_points[0]]
    intersection = test_segment_intersection(line, segment)

    if intersection:
        n_cuts += 1

    if n_cuts % 2 != 0 or is_inside is True or on_segment(point, segment):
        is_inside = True
    else:
        is_inside = False

    return is_inside, n_cuts

# SECTION 3
def section3_calculate_tangents(P, q):
    """
    Calculates the tangent points from a point 'q' to a convex polygon 'P'
    by detecting sign changes in the signed area.

    Args:
        P : list[list[float, float]]
            List of convex polygon vertices in order.
        q : list[float, float]
            Point from which tangents are to be calculated.

    Returns:
        list :
            List of polygon vertices that act as tangent points.
    """
    signs = []
    tangents = []

    for i in range(len(P) - 1):
        vi_1 = P[i]
        vi_2 = P[i+1]

        signs.append(sign(signed_area(vi_1, vi_2, q)))

        if len(signs) >= 2:
            if signs[i] != signs[i - 1]:
                tangents.append(vi_1)

    final_sign = sign(signed_area(vi_2, P[-1], q))
    if final_sign != signs[-1]:
        tangents.append(vi_2)

    return tangents

# SECTION 4
def polygon_area(p):
    """
    Calculates the signed area of a simple polygon given by its ordered vertices.
 
    Args:
        p : list[tuple[float, float]] : List of polygon vertices in counterclockwise order (or clockwise if negative area is desired).
 
    Returns:
        area : float : Signed area of the polygon. 
                       If vertices are in counterclockwise order, the area is positive.
                       If they are in clockwise order, it is negative.
    """
    area = 0
    for i in range(1, len(p)-1):
        area += signed_area(p[0], p[i], p[i+1])
    return area



# SECTION 5
def polygon_diagonals(P):
    """
    Calculates the internal and external diagonals of a simple polygon.
 
    Args:
        P : list[(float,float)] List of polygon vertices in counterclockwise order.
 
    Returns:
        (internal, external) : tuple[list, list] Two lists with internal and external diagonals, each formed by pairs of points ((x1,y1),(x2,y2)).
    """
 
    n = len(P)
    internal_diagonals = []
    external_diagonals = []
 
    for i in range(n):
        for j in range(i+1, n):
            # discard adjacent sides or polygon closure
            if abs(i-j) == 1 or (i == 0 and j == n-1):
                continue
 
            vi, vj = P[i], P[j]
            v_im1, v_ip1 = P[(i-1)%n], P[(i+1)%n]
            v_jm1, v_jp1 = P[(j-1)%n], P[(j+1)%n]
 
            # Check that it doesn't intersect any side
            intersects = False
            for k in range(n):
                a, b = P[k], P[(k+1)%n]
                if k in {i, j, (i-1)%n, (j-1)%n}:  # avoid adjacent sides
                    continue
                if test_segment_intersection((vi, vj), (a, b)):
                    intersects = True
                    break
            if intersects:
                continue  # not a valid diagonal, discard it
 
            # Local classification: convex or concave
            # vertex vi
            if signed_area(v_im1, vi, v_ip1) > 0:  # convex
                cond_i = (signed_area(vi, vj, v_im1) > 0) and (signed_area(vi, vj, v_ip1) < 0)
            else:  # concave
                cond_i = not ((signed_area(vi, vj, v_im1) < 0) and (signed_area(vi, vj, v_ip1) > 0))
            # vertex vj
            if signed_area(v_jm1, vj, v_jp1) > 0:  # convex
                cond_j = (signed_area(vj, vi, v_jm1) > 0) and (signed_area(vj, vi, v_jp1) < 0)
            else:  # concave
                cond_j = not ((signed_area(vj, vi, v_jm1) < 0) and (signed_area(vj, vi, v_jp1) > 0))
 
            # Final classification
            if cond_i and cond_j:
                internal_diagonals.append((vi, vj))
            else:
                external_diagonals.append((vi, vj))
 
    return internal_diagonals, external_diagonals


# SECTION 6
def concave_convex_vertices(p):
    """
    Classifies polygon vertices into convex and concave.

    Args:
        p : list[tuple[float, float]] :
            List of polygon vertices in order (clockwise or counterclockwise).

    Returns:
        convex : list[tuple[float, float]] :
            List of vertices that form turns with the same sign
            as the polygon orientation (convex vertices).
        
        concave : list[tuple[float, float]] :
            List of vertices that form turns with opposite sign
            to the polygon orientation (concave vertices).

    """
    orientation = polygon_orientation(p)
    convex, concave = [], []
    n = len(p)

    for i in range(n):
        a, b, c = p[i-1], p[i], p[(i+1) % n]
        turn = sign(signed_area(a, b, c))

        if turn == orientation:
            convex.append(b)
        elif turn == -orientation:
            concave.append(b)

    return convex, concave


# EXTRA SECTION - IN CASE OF NON-CONVEX POLYGON
def convex_hull_lee(points):
    """
    Calculates the convex hull of a set of points using the
    Lee algorithm (Gift Wrapping / Jarvis March).

    Args:
        points : list[tuple[float, float]]
            Set of points to wrap.

    Returns:
        list :
            Vertices of the convex hull in counterclockwise order.
            If there are fewer than 3 points, they are returned as is.
    """
    if len(points) < 3:
        return points

    start_point = y_min(points)

    hull = [start_point]
    current_point = start_point

    while True:
        next_point = points[0] if points[0] != current_point else points[1]

        for p in points:
            if p == current_point:
                continue

            pos = signed_area(current_point, next_point, p)

            if (pos > 0 or
                (pos == 0 and dist(current_point, p) > dist(current_point, next_point))):
                next_point = p

        current_point = next_point

        if current_point == start_point:
            break

        hull.append(current_point)


    return hull