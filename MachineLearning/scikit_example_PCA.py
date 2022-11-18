import numpy as np

from scipy import stats

e = np.exp(1)
np.random.seed(4)


def pdf(x):
    return 0.5 * (stats.norm(scale=0.25 / e).pdf(x) + stats.norm(scale=4 / e).pdf(x))


y = np.random.normal(scale=0.5, size=(30000))
x = np.random.normal(scale=0.5, size=(30000))
z = np.random.normal(scale=0.1, size=len(x))

density = pdf(x) * pdf(y)
pdf_z = pdf(5 * z)

density *= pdf_z

a = x + y
b = 2 * y
c = a - b + z

norm = np.sqrt(a.var() + b.var())
a /= norm
b /= norm



from sklearn.decomposition import PCA

import matplotlib.pyplot as plt

# unused but required import for doing 3d projections with matplotlib < 3.2
import mpl_toolkits.mplot3d  # noqa: F401


def plot_figs(fig_num, elev, azim):
    fig = plt.figure(fig_num, figsize=(4, 3))
    plt.clf()
    ax = fig.add_subplot(111, projection="3d", elev=elev, azim=azim)
    ax.set_position([0, 0, 0.95, 1])


    """
    a[::10] syntax:
    Simply means take 'a' starting from 0-th element and go the end skipping by 10.
    So, if 
    >>> ii=np.arange(100) 
    then
    >>> ii[::10]
    array([ 0, 10, 20, 30, 40, 50, 60, 70, 80, 90])
    The unmarked 1st and 2nd colons stand for 0-th and last elements
    
    Now this is generalized in this example
    >>> ii[7:33:3]
    array([ 7, 10, 13, 16, 19, 22, 25, 28, 31])
    i.e, start at element 7 of ii and go to element 33 skipping by 3
    Of course, here we stop at 31 since next element in sequnce would be 34
    """

    ax.scatter(a[::10], b[::10], c[::10], c=density[::10], marker="+", alpha=0.4)

    """
    np.c_[a,b,c] syntax:
    Here a,b,c are row vectors lenght N
    >>> Y=np.c_[a,b,c] 
    creates an array with N rows and 3 columns, where
    Y[irow,:] = [ a1.shaa[irow],b[irow],c[irow] ]
    """

    Y = np.c_[a, b, c]

    # Using SciPy's SVD, this would be:
    # _, pca_score, Vt = scipy.linalg.svd(Y, full_matrices=False)

    pca = PCA(n_components=3)
    pca.fit(Y)
    V = pca.components_.T

    x_pca_axis, y_pca_axis, z_pca_axis = 3 * V
    x_pca_plane = np.r_[x_pca_axis[:2], -x_pca_axis[1::-1]]
    y_pca_plane = np.r_[y_pca_axis[:2], -y_pca_axis[1::-1]]
    z_pca_plane = np.r_[z_pca_axis[:2], -z_pca_axis[1::-1]]
    x_pca_plane.shape = (2, 2)
    y_pca_plane.shape = (2, 2)
    z_pca_plane.shape = (2, 2)
    ax.plot_surface(x_pca_plane, y_pca_plane, z_pca_plane)
    ax.w_xaxis.set_ticklabels([])
    ax.w_yaxis.set_ticklabels([])
    ax.w_zaxis.set_ticklabels([])


elev = -40
azim = -80
plot_figs(1, elev, azim)

elev = 30
azim = 20
plot_figs(2, elev, azim)

plt.show()

