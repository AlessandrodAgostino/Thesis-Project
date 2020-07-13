import numpy as np
import itertools
import matplotlib.pyplot as plt
import os
import time

from scipy.spatial import Voronoi, ConvexHull, cKDTree
from SALib.sample import saltelli
from noise import pnoise2
from mpl_toolkits.mplot3d import Axes3D

from D3branch import *

def _inside_boundaries(vor, reg, boundaries):
    """
    This function check if a 'region' of a given 'Voronoi' tassellation lies
    inside some given 'boundaries'.

    Parameters:
        vor    : scipy.spatial.Voronoi object

        reg    : subset of vor.vertices to check

        boundaries : [ [min_x, max_x], [min_y, max_y], [min_z, max_z] ]

    Return: True if 'all' the verticies lies inside the region. False otherwise.
    """
    if -1 in reg: return False
    if len(reg)>1:
        reg = np.asarray([vor.vertices[pt] for pt in reg])
        if all(reg[:,0] > boundaries[0][0]) and all(reg[:,0] < boundaries[0][1]):
            if all(reg[:,1] > boundaries[1][0]) and all(reg[:,1] < boundaries[1][1]):
                if all(reg[:,2] > boundaries[2][0]) and all(reg[:,2] < boundaries[2][1]):
                    return True
        return False
    else: return False

def _plane_y_intersection(p1, p2, k=0):
    """
    Return the intersection between the plane of equation y = k and the line
    crossing 'p1' and 'p2'.
    """
    x = p1[0] + (k - p1[1]) / (p2[1] - p1[1]) * (p2[0] - p1[0])
    z = p1[2] + (k - p1[1]) / (p2[1] - p1[1]) * (p2[2] - p1[2])
    return np.asarray((x,k,z))

def _box_and_spheres(ramification, y = 0):
    """
    Given a ramification this function returns a list spatial information.

    Parameters:
        ramification: The ramification to work onto.

        y:            Section plane's height.

    Returns a tuple with:
        boundaries: The box interval that contains all the ramification, with an
        empty padding proportional to free ends' spheres' radius.

        small_boundaries: A smaller volume adjacent to section plane.
        (MAYBE INSERT THICKNESS F THE VOLUME)

        int_spheres: The list of the spheres that lie near the section plane.
        Those are the interesting spheres.

        sph_rad: The radius of the spheres.
    """
    #Extracting free end's spheres and radius
    spheres  = []
    sph_rad  = 0
    max_iter = np.log2((len(ramification)+1)) - 1

    for br in ramification:
        if br.iter_lev == max_iter:
            cent = br.pos + vector(*(br.length * br.drct))
            spheres.append(np.asarray((cent.x, cent.y, cent.z)))
            if not sph_rad: sph_rad = br.length

    #Interesting spheres. Those around x,y-plane
    int_spheres = [sph for sph in spheres if (sph[1] >  y - 2*sph_rad) & (sph[1] <  y + 2*sph_rad)]

    #Boundaries for random sampling in a volume with a padding proportional to spheres' radius
    max_box = np.max(np.asarray(spheres), axis = 0) + sph_rad*2
    min_box = np.min(np.asarray(spheres), axis = 0) - sph_rad*2

    boundaries = np.array([[min_box[0], max_box[0]],
                          [min_box[1], max_box[1]],
                          [min_box[2], max_box[2]]])

    small_boundaries = np.array([[min_box[0], max_box[0]],
                                [y - sph_rad,  y + sph_rad],
                                [min_box[2], max_box[2]]])

    return boundaries, small_boundaries, int_spheres, sph_rad

def _get_nuclei_radius(boundaries, N_points, nuclei_scale = 20):
    """
    Returns the estimate nuclei's radius.
    """
    box_volume = np.prod([ max - min for min, max in boundaries])
    cell_dim = np.cbrt(box_volume/N_points)
    nuclei_rad = cell_dim/nuclei_scale
    return nuclei_rad

def _get_vor_points(boundaries,  N_points, small_boundaries = None, method = 'saltelli'):
    """
    This function defines the problem for a low discrepancy random sampling
    inside the volume defined by small_boundaries.

    Parameters:
        boundaries: Volume into which make the external sampling

        small_boundaries: Volume into which select points

        N_points: Number of points to be filtered inside small_boundaries

    Return:
        points: All the filtered points as np.array
    """
        #Mapping the points to the Boundaries
    if method == 'saltelli':
        if small_boundaries == None:
            small_boundaries = boundaries

        problem = {'num_vars': 3,
                   'names': ['x', 'y', 'z'],
                   'bounds': boundaries}
        N = 10
        points = []
        while len(points) < N_points:
            N = N*2
            points = saltelli.sample(problem, N)
            points = points[(points[:,1] > small_boundaries[1][0]) & (points[:,1] < small_boundaries[1][1])]
            points = points[:N_points]

    elif method == 'r-sequence':
        #Generating points through a recurrence-rule
        d = 3
        s_0 = 0.5
        x=2.0000
        for i in range(10):
            x = pow(1+x,1/(d+1))
        g = x

        alphas = np.power(1/g, np.arange(1,d+1))
        points = np.mod(np.outer(alphas,np.arange(1, N_points+1)) + s_0, 1)

        #Mapping the points to the Boundaries
        lengths = small_boundaries[:,1] - small_boundaries[:,0]
        points = ((points.T * lengths) + small_boundaries[:,0])

    elif method == 'lattice':
        n_sample = np.floor(np.power(N_points, 1/3)).astype(int)
        coords = np.zeros((3,n_sample))
        coords[0] = np.linspace(small_boundaries[0][0], small_boundaries[0][1], n_sample)
        coords[1] = np.linspace(small_boundaries[1][0], small_boundaries[1][1], n_sample)
        coords[2] = np.linspace(small_boundaries[2][0], small_boundaries[2][1], n_sample)
        points = np.asarray([ pt for pt in itertools.product(coords[0], coords[1], coords[2])])

    elif method == 'uniform':
        points = np.random.uniform(size = N_points*3)
        points = points.reshape(-1,3)

        lengths = small_boundaries[:,1] - small_boundaries[:,0]
        points = ((points * lengths) + small_boundaries[:,0])

    else: print('Unknow method')

    return points

def _draw_section(vor, cropped_reg, region_id, palette, h, nuclei_rad, draw_nuclei):
    """
    This auxiliary funtion creates the plt.figure of a section of a certain Voronoi
    tassellation at an height of h.

    Parametrs:
        vor: Voronoi tasselation to section off.

        cropped_reg: regions of interest

        region_id: identity of the interesting regions

        palette: color map for the identity

        h: height of the section

    Returns:
        plt.figure representing the section
    """
    fig = plt.figure(figsize=(8,8))
    plt.gca().get_yaxis().set_ticks([])
    plt.gca().get_xaxis().set_ticks([])
    plt.gca().set_facecolor("grey")

    #Nuclei Projection
    if draw_nuclei:
        projectable_nuclei = [pt for pt in vor.points if (np.abs(pt[1] - h) <  nuclei_rad*10) ]
        circles = [(pt[0], pt[2], nuclei_rad) for pt in projectable_nuclei]

        # projectable_nuclei = [pt for pt in vor.points if (np.abs(pt[1] - h) <  nuclei_rad) ]
        # circles = [(pt[0], pt[2], np.sqrt(nuclei_rad**2 - (np.abs(pt[1] - h ))**2)) for pt in projectable_nuclei]

        for circ in circles:
            circle = plt.Circle((circ[0], circ[1]), circ[2], color = palette[3], alpha = 0.7)
            plt.gca().add_artist(circle)

    #Draw Cells
    for n, reg in enumerate(cropped_reg):
        ind_abo = [ ver for ver in vor.vertices[reg] if ver[1] > h]
        ind_bel = [ ver for ver in vor.vertices[reg] if ver[1] <= h]

        couples = list(itertools.product(ind_abo,ind_bel))
        if couples:
            intersection_point = [ _plane_y_intersection(v1, v2, k = h) for v1, v2 in couples]
            intersection_point = np.asarray(intersection_point)
            drawing_hull = ConvexHull(intersection_point[:,[0,2]]).vertices
            plt.gca().fill(intersection_point[drawing_hull,0],
                    intersection_point[drawing_hull,2],
                    color = palette[region_id[n]])
    return fig

def _draw_section_2(vor, cropped_reg, palette, h, nuclei_rad, draw_nuclei, bound):
    """
    This auxiliary funtion creates the plt.figure of a section of a certain Voronoi
    tassellation at an height of h.

    Parametrs:
        vor: Voronoi tasselation to section off.

        cropped_reg: regions of interest

        region_id: identity of the interesting regions

        palette: color map for the identity

        h: height of the section

    Returns:
        plt.figure representing the section
    """
    fig = plt.figure(figsize=(8,8))
    plt.gca().get_yaxis().set_ticks([])
    plt.gca().get_xaxis().set_ticks([])
    plt.gca().set_xlim(*bound[0])
    plt.gca().set_ylim(*bound[2])
    plt.gca().set_facecolor("grey")

    circles = []

    patches = []
    for reg in cropped_reg:
        ind_abo = [ ver for ver in reg['vertices'] if ver[1] > h]
        ind_bel = [ ver for ver in reg['vertices'] if ver[1] <= h]

        couples = list(itertools.product(ind_abo,ind_bel))
        if couples:
            intersection_point = [ _plane_y_intersection(v1, v2, k = h) for v1, v2 in couples]
            intersection_point = np.asarray(intersection_point)
            drawing_hull = ConvexHull(intersection_point[:,[0,2]]).vertices

            hull_points = np.stack((intersection_point[drawing_hull,0],intersection_point[drawing_hull,2]), axis = -1)
            from matplotlib.patches import Polygon
            from matplotlib.collections import PatchCollection

            if reg['identity'] == 2:
                polygon = Polygon(hull_points, True, ec = palette[2], fc = (1,1,1,1), lw =1)

            else:
                polygon = Polygon(hull_points, True, fc = palette[reg['identity']])
            plt.gca().add_artist(polygon)

            #Nuclei
            if reg['identity'] < 2:
                circles.append((reg['nucleus'][0], reg['nucleus'][2], nuclei_rad))

    if draw_nuclei:
        for circ in circles:
            circle = plt.Circle((circ[0], circ[1]), circ[2], color = palette[2], alpha = 0.7)
            plt.gca().add_artist(circle)

    return fig

def _draw_noise(RGB = True, cmap = 'Purples', noise_density = 20):
    """
    This funtion produces a random Perlin noise image

    Parameters:
        cmap: color map for noise values
        noise_density: determines the fineness of the noise

    Return:
        plt.figure objet of noise
    """
    #Settings
    density = 1000
    offset = np.random.randint(1000, size =1)

    vec_noise = np.vectorize(pnoise2)

    fig = plt.figure(figsize=(8,8))
    plt.gca().get_yaxis().set_ticks([])
    plt.gca().get_xaxis().set_ticks([])
    plt.gca().set_facecolor("grey")
    plt.tight_layout()

    if RGB:
        x_val  = np.linspace(offset, 3*noise_density + offset, 3*density)
        y_val  = np.linspace(0, noise_density, density)
        x, y = np.meshgrid(x_val, y_val, sparse = False)
        z = vec_noise(x, y) * 1/6 + 1/2
        z_stack = np.dstack([z[:,:density], z[:,density:2*density], z[:,2*density:]])
        plt.imshow(z_stack)

    else:
        #Creating Perlin noise
        val  = np.linspace(offset, noise_density + offset, density)
        x, y = np.meshgrid(val, val, sparse = False)
        z = vec_noise(x, y) * 128 + 128
        #Drawing contour
        plt.contourf(x, y, z, cmap = cmap)

    return fig

def _create_surf():
    x_max = 1.5
    y_max = 1.5
    n_pts = 2000
    scale = 5 #Delicate Parameter
    offset = 10 #Acts like a seed

    ED_height = 1
    EK_height = 0.2

    K_thick = 0.4
    E_thick = 0.5
    D_thick = 0.5

    vec_noise = np.vectorize(pnoise2)

    x_val  = np.linspace(offset, x_max + offset, n_pts)
    y_val  = np.linspace(offset, y_max + offset, n_pts)
    xx, yy = np.meshgrid(x_val, y_val, sparse = False)
    zz = vec_noise(xx, yy) * EK_height + E_thick

    print(f"EK surf in range [{zz.min():.2f},{zz.max():.2f}]")

    xx = (xx - offset) * scale + offset
    yy = (yy - offset) * scale + offset
    zz1 = vec_noise(xx, yy) * ED_height
    print(f"ED surf in range [{zz1.min():.2f},{zz1.max():.2f}]")

    x_grain = np.unique(xx)[1] - np.unique(xx)[0]
    y_grain = np.unique(yy)[1] - np.unique(yy)[0]

    disc_grain = (x_grain, y_grain)

    EK_bound = zz
    ED_bound = zz1
    K_bound = EK_bound + K_thick

    max_z = K_bound.max() + 2*K_thick
    min_z = -ED_height - D_thick

    boundaries = np.array([[0     , x_max * scale ],
                           [0     , y_max * scale ],
                           [min_z , max_z ]])


    # fig = plt.figure(figsize=(15,10))
    # ax = fig.gca(projection='3d')
    # surf = ax.plot_trisurf(*(cc.flatten() for cc in ED_bound), label='Epid - Derm')
    # surf1 = ax.plot_trisurf(*(cc.flatten() for cc in EK_bound), label='Kera - Epid')
    # surf._facecolors2d=surf._facecolors3d
    # surf._edgecolors2d=surf._edgecolors3d
    # surf1._facecolors2d=surf1._facecolors3d
    # surf1._edgecolors2d=surf1._edgecolors3d
    # leg = ax.legend()
    # ax.set_title("Perlin Noise", fontsize=20)

    return boundaries, ED_bound, EK_bound, K_bound, disc_grain

def _surf_identity(nuclei, bound):
    identity  = []
    for xn, yn, zn in nuclei:
        x = xn.astype('int')
        y = yn.astype('int')
        identity.append(zn > bound[x,y]) #True: above, False: beneath
    return identity

def derma_section(N_points=600, n_slices = 1, saving_path='.', dpi = 100, draw_nuclei=True):
    plane_distance = 0 #????
    boundaries, ED_bound, EK_bound, K_bound, (x_grain, y_grain)= _create_surf()
    vor_points = _get_vor_points(boundaries, N_points) #Getting point for creating Voronoi tassellation
    vor = Voronoi(vor_points) #Creating the tassellation

    nuclei_rad = _get_nuclei_radius(boundaries, N_points)
    nuclei = np.copy(vor.points)
    nuclei[:,0:2] = np.floor_divide(nuclei[:,0:2], np.array(x_grain, y_grain))

    K_identity  = _surf_identity(nuclei, K_bound)  # Void
    EK_identity = _surf_identity(nuclei, EK_bound) # Void - Keratine
    ED_identity = _surf_identity(nuclei, ED_bound) # Void - Keratine - Epidermis

    nuclei_id = np.asarray(ED_identity)*1 + np.asarray(EK_identity)*1 + np.asarray(K_identity)*1
    # 0:Dermis, 1:Epidermis, 2:Keratine, 3:Void

    cropped_reg = []
    for n,pt in enumerate(vor.point_region):
        reg = vor.regions[pt]
        if _inside_boundaries(vor, reg, boundaries):
            cr_reg = {'nucleus' : vor.points[n],
                      'vertices': [vor.vertices[id] for id in reg],
                      'identity': nuclei_id[n]}  # 0:Dermis, 1:Epidermis, 2:Keratine, 3:Void
            cropped_reg.append(cr_reg)

    l = np.arange(0,n_slices)
    slices = (l - np.floor(len(l)/2)).astype('int8')

    #Plotting settings
    turquoise = color.hsv_to_rgb(vector(0.5,1,0.8))
    red       = color.red #Some colors
    white     = color.white
    black     = color.black
    colors    = [red, turquoise, white,  black] # 0:Dermis, 1:Epidermis, 2:Keratine, 3:Void

    IDENTITIES = {'Dermis':'red', 'Epidermis':'turquoise', 'Keratine':'white', 'Void':'black'}
    for k, v in IDENTITIES.items():
        print(f'{k[:5]}: \t {v}')

    Figures   =  [] #List to which append all the drawings

    # 3D Plotting
    # scene     = canvas(width=800, height=600, center=vector(5,5,0), background=color.gray(0.6))
    # for reg in cropped_reg[:1000]:
    #     conv_hull= ConvexHull(reg['vertices'])
    #     simpl = []
    #     for sim in conv_hull.simplices:
    #         pts = [conv_hull.points[pt] for pt in sim]
    #         simpl.append( triangle( vs=[vertex( pos     = vector(*ver),
    #                                             color   = colors[reg['identity']],
    #                                             opacity = 0.2) for ver in pts]))

    #Loop for Drawing and Saving every SLICE
    #Previous palette
    palette = [[0.81176471, 0.62745098, 0.78039216],
               [0.49803922, 0.34509804, 0.58823529],
               [0.34,       0.19,       0.63],
               [0.9254902,  0.89411765, 0.91372549]]
    y = 0.5
    for n_s in slices:
        dy = plane_distance * n_s
        fig = _draw_section_2(vor, cropped_reg, palette, draw_nuclei = draw_nuclei, h = y+dy, nuclei_rad = nuclei_rad, bound = boundaries)
        fig.savefig(os.path.join(saving_path + f'N_{N_points}_slice.png'),
                    #bbox_inches='tight',
                    dpi=dpi)
        plt.close(fig)

    # #Drawing the LABEL image
    # fig = _draw_section(vor, cropped_reg, region_id, lab_colors, draw_nuclei = draw_nuclei, h = y, nuclei_rad = nuclei_rad)
    # fig.savefig(os.path.join(saving_path + f'N_{N_points}_slice.png'),
    #             #bbox_inches='tight',
    #             dpi=dpi)
    plt.close(fig)


def pancr_section( iteration_level = 3,rotation = False, seed = None, y = 0, N_points = 5000, n_slices = 3, saving_path = '', noise_density = 20, plane_distance = 0.05, draw_nuclei = True, sampling_method = 'saltelli'):
    """
    This function draws and saves the images resulting from the slicing of a ramification.

    Parameters:
        iteration_level: The # of iterative biforcation made in the ramification

        rotation: If the Tree has to be rotated or not in a random direction

        seed: The seed for the random direction

        y: The height on the y axis for the section

        N_points: # of points in the volme adjacent to the section plane, hence the density on points. N_points > 5000 is suggested

        n_slices: # of slices to be made around the section plane

        saving_path: Where to store images

        noise_density: Parameter that regulates Perlin noise density

        plane_distance: Distance between different section planes


    Returns:
        times: The complete time report for the section process.
    """
    Ramification = createTree(iter = iteration_level, rotation = rotation, seed = seed) #Creating the ramification object
    boundaries, small_boundaries, int_spheres, sph_rad = _box_and_spheres(Ramification, y) #Getting spatial informations
    vor_points = _get_vor_points(boundaries, small_boundaries, N_points, method = sampling_method) #Getting point for creating Voronoi tassellation

    vor = Voronoi(vor_points) #Creating the tassellation
    cropped_reg = [ reg for reg in vor.regions if _inside_boundaries(vor, reg, small_boundaries)] #Cropping out the regions that lies outside the boundaries
    nuclei_rad = _get_nuclei_radius(boundaries, N_points)

    tree = cKDTree(vor.vertices) #Creating a Tree object for fast distances computation
    #Vertices identity
    vertices_truth = np.zeros(len(vor.vertices)) #Array where to store identities
    try:
        inside_indexes = tree.query_ball_point(int_spheres, sph_rad)
        for ind in inside_indexes:
            vertices_truth[ind] = 1
        Null_Image = False
    except:
        Null_Image = True
        print(f'Seed {seed} gave uninterseting section')
    vertices_truth = vertices_truth.astype(int) #0: Outside, 1: Inside

    #Region identity
    region_id = np.zeros(len(cropped_reg))
    for n, reg in enumerate(cropped_reg):
        region_id[n] = any(vertices_truth[reg]) + all(vertices_truth[reg])
    region_id = region_id.astype(int) # 0:Outside, 1:Partially, 2:Inside

    #DRAWING and SAVING SETINGS
    #--------------------------------------------------------------------------
    cmap = plt.get_cmap('Purples')
    #palette = [cmap(0.1)[:-1], cmap(0.8)[:-1], cmap(0.4)[:-1]]
    lab_colors = ['w', 'c', 'r', 'y']
    dpi = 100
    #Previous palette
    palette = [[0.9254902,  0.89411765, 0.91372549],
               [0.49803922, 0.34509804, 0.58823529],
               [0.81176471, 0.62745098, 0.78039216],
               [0.34,       0.19,       0.63]]

    l = np.arange(0,n_slices)
    slices = (l - np.floor(len(l)/2)).astype('int8')
    #Loop for Drawing and Saving every SLICE
    for n_s in slices:
        dy = plane_distance * n_s
        fig = _draw_section(vor, cropped_reg, region_id, palette, draw_nuclei = draw_nuclei, h = y+dy, nuclei_rad = nuclei_rad)
        fig.savefig(os.path.join(saving_path + f'N_{N_points}_seed_{seed}_sl_{n_s}.png'),
                    #bbox_inches='tight',
                    dpi=dpi)
        plt.close(fig)

    #Drawing the LABEL image
    fig = _draw_section(vor, cropped_reg, region_id, lab_colors, draw_nuclei = draw_nuclei, h = y, nuclei_rad = nuclei_rad)
    fig.savefig(os.path.join(saving_path + f'N_{N_points}_seed_{seed}_label_uncropped.png'),
                #bbox_inches='tight',
                dpi=dpi)
    plt.close(fig)

#%%
def different_density_benchmarks():
    MAX = 45000
    MIN = 35000
    STEP = 5000
    copies = 10
    n_slices = 4

    seeds = (s for s in np.random.randint(1000, size= int((MAX - MIN) / STEP * copies *2)))
    time_df = pd.DataFrame()

    for N_points in np.arange(MIN, MAX, STEP):
        for n in range(copies):
            fig, times = section(iteration_level = 3,
                                 y = -2,
                                 n_slices = n_slices,
                                 rotation = True,
                                 seed = next(seeds),
                                 N_points = N_points,
                                 saving_path = '')

            time_df = pd.concat([time_df, pd.DataFrame(times)], ignore_index = True)
            time_df.to_csv('Times/time_measures.csv')
        print(f'{N_points} figures written.')

if __name__ == '__main__':
    derma_section(N_points = 120000)
