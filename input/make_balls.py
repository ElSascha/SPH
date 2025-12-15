'''Makes two spheres of constant particle density that slowly move to each other.
See constants lower down to set specifics.'''
import numpy as np
import matplotlib.pyplot as plt


def points_on_sphere_surface(N: int):
    # https://stackoverflow.com/a/26127012
    # from https://arxiv.org/pdf/0912.4540
    # adapted for numpy
    '''Isotropically distributes `N` points on the surface of a unit sphere
    centered on the origin, using the golden ratio.'''
    # golden ratio
    phi = np.pi * (np.sqrt(5.) - 1.)
    # equally spaced along y
    i = np.arange(N)
    y = 1 - (i / float(N-1)) * 2
    y = np.clip(y, -1, 1)
    # derive x and z from y and angle theta (in steps of golden ratio)
    r = np.sqrt(1 - y * y)
    theta = phi * i
    x = np.cos(theta) * r
    z = np.sin(theta) * r
    return np.stack([x, y, z])


def total_points_in_ball(radius: float, target_density: float) -> int:
    '''Returns the number of particles needed to fill a ball of radius `radius`
     at average density `target_density`. '''
    volume = 4 * np.pi / 3 * radius**3
    return np.round(volume * target_density, 0).astype(int)


def points_union(point_arrays: list[np.ndarray]) -> np.ndarray:
    '''Concatenates the arrays along the rows (not the columns).'''
    return np.concatenate(point_arrays, axis=1)


def make_ball(radius: float, delta_r: float):
    '''Creates a solid ball of particles, each roughly `delta_r` away from its
    neighbours, with total radius `radius`, centered on the origin.
    
    Expect number of particles ~= 4/3 pi (radius/delta_r)^3'''
    # evaluating total points per shell
    N_total = total_points_in_ball(radius, delta_r**-3)
    shells = np.arange(delta_r, radius+delta_r, delta_r)
    shell_distribution = shells**2
    shell_distribution *= N_total / np.sum(shell_distribution)
    # making shells
    shell_points = [
        points_on_sphere_surface(N) * shell
        for shell, N in zip(shells, shell_distribution)
    ]
    # full ball
    return points_union(shell_points)


def make_ball_cubic(radius: float, delta_r: float):
    '''Makes a ball of radius `radius` at the origin with points spaced roughly
    `delta_r` apart in a simple cubic grid.'''
    spacing = delta_r # can be more complex geometrically for consistent density, but doing it simple for now
    grid_points = np.mgrid[-radius:radius:spacing,-radius:radius:spacing,-radius:radius:spacing]
    grid_points = grid_points.reshape(3, -1)
    ball_points = grid_points[:, np.linalg.norm(grid_points, axis=0) < radius]
    return ball_points


def make_ball_hcp(radius: float, delta_r: float):
    '''Makes a ball of radius `radius` at the origin with points spaced roughly
    `delta_r` apart in a hexagonal close packed grid.'''
    # math: https://en.wikipedia.org/wiki/Close-packing_of_equal_spheres
    r = delta_r/2
    N_1d = np.ceil(radius*2/delta_r).astype(int)
    hex_row_sep = np.sqrt(3)
    hex_plane_sep = 2*np.sqrt(6)/3
    index_grid = (np.mgrid[0:N_1d, 0:N_1d, 0:N_1d]).reshape(3, -1)
    i, j, k = index_grid
    x = 2*i + np.mod(j + k,2)
    y = hex_row_sep * (j + np.mod(k, 2) / 3)
    z = hex_plane_sep * k
    points = np.stack([x * r, y * r, z * r])
    # shift to center on origin
    points -= np.mean(points, axis=1)[:, None]
    # cull points outside ball
    ball = points[:, np.linalg.norm(points, axis=0) < radius]
    return ball


def shift_points(points: np.ndarray, dx=0.0, dy=0.0, dz=0.0):
    '''Move all points in `points` by the specified x, y, and z shifts.'''
    shift = np.asarray([dx, dy, dz], float)
    return points[:, :] + shift[:, None]


def show(a: np.ndarray):
    '''Plots and displays a figure of the points in a.'''
    fig = plt.figure()
    ax: plt.Axes = fig.add_subplot(projection='3d')
    ax.scatter(a[0, :], a[1, :], a[2, :], s=10.0)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.tight_layout()
    ax.set_aspect('equal')
    plt.show()


def random_from_cdf(cdf: tuple[np.ndarray, np.ndarray], N: int):
    '''Draws `N` random values from the distribution specified by `cdf`, whose
    first array is the domain to draw from, and the second array is the
    cumulative distribution function.'''
    x, y = cdf
    # draw uniform on domain [0, 1]
    u = np.random.random(N)
    # find cdf indices
    i0 = np.clip(np.searchsorted(y, u, side='left')-1, 0, len(y)-1)
    i1 = np.clip(np.searchsorted(y, u, side='right'), 0, len(y)-1)
    # lerp to x domain
    x0, x1, y0, y1 = x[i0], x[i1], y[i0], y[i1]
    return x0 + (u-y0) * (x1-x0)/(y1-y0)


def make_impact_flaw_cdf():
    '''Makes a `cdf` for the distribution seen in the flaw activation energies
    of the impact example in miluphcuda.'''
    xlow, xhigh = 1.4e-4, 2.1e-4
    x = np.linspace(xlow, xhigh, 20)
    pdf = (x-xlow)**2
    pdf /= np.sum(pdf)
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]
    return (x, cdf)


impact_flaw_cdf = make_impact_flaw_cdf()


def points_to_particles(
    points: np.ndarray,
    velocity: tuple[float, float, float],
    mass: float,
    density: float,
    h: float = 0.5,
    distention: float = 1.25,
    flaw_activation_energy_dist_cdf: tuple[np.ndarray, np.ndarray] = impact_flaw_cdf
):
    '''Assigns the particles at points `points` the provided properties. All
    other properties (not counting default arguments of this function) are
    defaulted to 0.0.'''
    N = points.shape[-1]
    p = np.zeros((24, N))
    p[:3, :] = points
    p[3:6, :] = np.repeat(velocity, N).reshape(3, N)
    p[6, :] = mass
    p[7, :] = density
    p[9, :] = h
    p[22, :] = distention
    return p


# units are SI unless otherwise specified
AGGREGATE_RADIUS = 0.1 # 10 cm
AGGREGATE_SHIFT = AGGREGATE_RADIUS * 2 # 40% of R distance between balls
POINT_DISTANCE = AGGREGATE_RADIUS / 4 # effect: 4 shells (but still ~isotropic distribution)

NR_INTERACTIONS_NORMAL = 50 # in the middle of the aggregate, expect N interaction partners (10% of max=512)
SMOOTHING_LENGTH = POINT_DISTANCE * (NR_INTERACTIONS_NORMAL)**(1/3)

DISTENTION = 1.00
RHO_0 = 2700.0 # monomer material density
PARTICLE_DENSITY = RHO_0 / DISTENTION
PARTICLE_MASS = PARTICLE_DENSITY * POINT_DISTANCE**3 # mass/particle = mass density / number density

v = 1.0 # 100 cm/s = 1/ms

#print(AGGREGATE_RADIUS, POINT_DISTANCE, SMOOTHING_LENGTH, PARTICLE_DENSITY, PARTICLE_MASS)
#make_ball_hcp(AGGREGATE_RADIUS, POINT_DISTANCE)
# test
if False:
    cubic = make_ball_cubic(AGGREGATE_RADIUS, POINT_DISTANCE)
    hcp = make_ball_cubic(AGGREGATE_RADIUS, POINT_DISTANCE)
    fib = make_ball(AGGREGATE_RADIUS, POINT_DISTANCE)

    cubic = shift_points(cubic, -AGGREGATE_SHIFT*5)
    hcp = shift_points(hcp, AGGREGATE_SHIFT*5)
    show(points_union([hcp, fib, cubic]))

# make geometry
left = make_ball(AGGREGATE_RADIUS, POINT_DISTANCE)
right = left.copy()
left = shift_points(left, dx=-AGGREGATE_SHIFT)
right = shift_points(right, dx=AGGREGATE_SHIFT)
print(left.shape[-1]*2)
show(points_union([left, right]))

# set particle properties
p_left = points_to_particles(
    left, (v, 0, 0), PARTICLE_MASS, PARTICLE_DENSITY, h=SMOOTHING_LENGTH, distention=DISTENTION
)
p_right = points_to_particles(
    right, (-v, 0, 0), PARTICLE_MASS, PARTICLE_DENSITY, h=SMOOTHING_LENGTH, distention=DISTENTION
)
p = points_union([p_left, p_right])

# output in HDF5 format matching create_distribution.py
import h5py

# Extract positions, velocities, mass, smoothing_length
positions = p[:3, :].T  # shape (N, 3)
velocities = p[3:6, :].T  # shape (N, 3)
mass = p[6, :]  # shape (N,)
smoothing_length = p[9, :]  # shape (N,)

with h5py.File('stable_ball.h5', 'w') as f:
    f.create_dataset('positions', data=positions)
    f.create_dataset('velocities', data=velocities)
    f.create_dataset('mass', data=mass)
    f.create_dataset('smoothing_length', data=smoothing_length)
print(f"Created HDF5 file with {positions.shape[0]} particles: 'stable_ball.h5'")
