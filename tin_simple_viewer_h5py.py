
#from rasputin import tin_repository as tin
# instead of rasputin we use pip available h5py to parse h5 files
import h5py
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib import colors
# from matplotlib import cm
import numpy as np
import matplotlib.tri as mtri
import pandas as pd
from pyproj import Proj
import seaborn as sns
import matplotlib.patches as mpatches
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import scipy
from scipy import stats
from sklearn.preprocessing import minmax_scale


def get_land_type(rgb, land_covers_available)->(int,str):
    for (v,n,r,g,b) in land_covers_available:
        rgbs = (r,g, b)
        (value,name) = (v,n)
        if rgb == rgbs:
            return (value,name, rgbs)
        else:
            pass

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]

def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list.

        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.

        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0, 1, len(rgb_list)))

    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = colors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp

class TinSimpleViewerH5PY:
    """
    Simple Tin Viewer
    Uses H5PY to parse the Tin-mesh file created by rasptuin software (https://github.com/expertanalytics/rasputin) ,
    converts to latlon and plots with matplotlib tools
    """

    def __init__(self):
        self._arr = []
        self._effarr = []
        self._elevarr = []
        self._latarr = []
        self._lonarr = []
        self._zfaces = []
        self._pv = []
        self._fv = []
        self._cv =[]
        self._cv_norm = []
        self._materials = []

    def lonlat(self,xarr, yarr,projection="+proj=utm +zone=33, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs"):
        latlon_m = {'lat': xarr,
                    'lon': yarr}
        #myProj = Proj("+proj=utm +zone=33, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
        myProj = Proj(projection)
        lon, lat = myProj(latlon_m['lat'], latlon_m['lon'], inverse=True)
        return lon, lat

    def distance(self,x1, y1, z1, x2, y2, z2):
        return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)

    def areaf(self,x1, y1, z1, x2, y2, z2, x3, y3, z3):
        a = self.distance(x1, y1, z1, x2, y2, z2)
        b = self.distance(x2, y2, z2, x3, y3, z3)
        c = self.distance(x3, y3, z3, x1, y1, z1)
        s = (a + b + c) * 0.5
        return np.sqrt(s * (s - a) * (s - b) * (s - c))


    def parse_tins(self,tin_data_folder,tid):

        pathtofile = Path(tin_data_folder+"/"+tid+".h5")

        #tin_repo = tin.TinRepository(path=Path(tin_data_folder))
        tin_repo = h5py.File(pathtofile, "r")
        # from here I already know teh structure of the file, by asking tin_repo.keys() one can check the main groups
        # in the rasputin generated tins main groups are "information" and "tin"
        tins = tin_repo.get("tin")
        # tins.keys() return "face_fields", "faces", "points"
        points = tins.get("points") # this is already a dataset, one cah check points.shape,
        # in the rasputin generated tin the shape is normally (L,3), where L -- is number of vertexes
        xarr = points[:,0]
        yarr = points[:,1]
        zarr = points[:,2]
        faces = tins.get('faces')
        # now we got all our vertexes, but still need to gain information about land_types
        # colors and codes for land_cover type are kept with tins
        face_fields = tins.get("face_fields") # this one carries additional information about cover_color and cover_type
        tincolors = face_fields.get("cover_color") #shape is same as for faces, as it is a color of face
        land_cover_type = face_fields.get("cover_type")
        # but text information about the name of the land_type is kept in information group
        info = tin_repo.get("information")
        # info.keys() returns "land_covers"
        land_covers = info.get("land_covers")

        materials = [(int(np.round(r * 255)), int(np.round(g * 255)), int(np.round(b * 255))) for (r, g, b) in
                     tincolors]
        pv = np.asarray(points)
        fv = np.asarray(faces)
        cv = np.asarray(tincolors)
        cv_norm = []
        for (r,g,b) in cv:
            cv_norm.append((r/255,g/255,b/255))

        # projstring = info["tin"]["projection"] # this is not parsed by h5py, not sure why
        projstring = "EPSG:32645"  # default to Nepal
        # print(materials)
        # print(land_cover_type)
        # print(land_covers)
        lt_values = []
        lt_colors = []
        for c in materials:
            (value, name, ltc) = get_land_type(c, land_covers)
            lt_values.append(value)
            lt_colors.append(ltc)
            # print(name)
        return pv, fv, cv, cv_norm, materials, lt_values, lt_colors, land_covers, projstring, land_cover_type

    def view_tins_terrain(self,tin_data_folder, tid, jpg_folder):
        # parse tins
        self._pv, self._fv, self._cv, self._cv_norm,self._materials, lt_values, lt_colors, land_covers, projstring, land_cover_type = self.parse_tins(tin_data_folder, tid)
        # set up figure
        sns.set(font_scale=1.6)
        fig = plt.figure(figsize=(12, 7))
        ax = fig.add_subplot(111)
        numtins = 0
        patches = []
        lonlistarr = []
        latlistarr = []
        elevlistarr = []
        areaslistarr = []
        effareaslistarr = []
        xarr = []
        yarr = []
        zarr = []
        for c in self._pv:
            xarr.append(c[0])
            yarr.append(c[1])
            zarr.append(c[2])
        self._lon, self._lat = self.lonlat(xarr, yarr,projstring)
        lonlistarr.append(xarr)
        latlistarr.append(yarr)
        elevlistarr.append(zarr)
        area = []
        effarea = []
        zfarr = []
        for k in self._fv:
            v0 = self._pv[k[0]]
            v1 = self._pv[k[1]]
            v2 = self._pv[k[2]]
            self._arr.append(self.areaf(v0[0], v0[1], v0[2], v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]))
            self._effarr.append(self.areaf(v0[0], v0[1], 0.0, v1[0], v1[1], 0.0, v2[0], v2[1], 0.0))
            mid_point = np.asarray(
                    [(v0[0] + v1[0] + v2[0]) / 3, (v0[1] + v1[1] + v2[1]) / 3, (v0[2] + v1[2] + v2[2]) / 3])
            zfarr.append(mid_point[2])
        self._zfaces = np.asarray(zfarr)

        areaslistarr.append(np.asarray(area))
        effareaslistarr.append(np.asarray(effarea))

        plt.tripcolor(self._lon, self._lat, self._fv, facecolors=self._zfaces, edgecolors='k', cmap="terrain") # plot terrain
        # plt.tripcolor(lon, lat, fv,  facecolors=self._zfaces, edgecolors='k', cmap='terrain', vmin = 0, vmax = 580)

        print(len(self._zfaces))
        print(len(self._fv))

        numtins += len(self._fv)
        title = 'Ncells: ' + str(numtins)
        plt.xlabel('Longitude, [deg]', fontsize=16)
        plt.ylabel('Latitude, [deg]', fontsize=16)
        plt.title(title, fontsize=18)
        plt.colorbar()
        print(land_covers)
        self._elevarr = np.array([elem for singleList in elevlistarr for elem in singleList])
        figname = jpg_folder + '/' + tid + '.jpg'
        #fig.savefig(figname)
        plt.show()


    def view_tins_land_cover(self,tin_data_folder, tid, jpg_folder):
        # parse tins
        self._pv, self._fv, self._cv, self._cv_norm,self._materials, lt_values, lt_colors, land_covers, projstring, land_cover_type = self.parse_tins(tin_data_folder, tid)
        # set up figure
        sns.set(font_scale=1.6)
        fig = plt.figure(figsize=(12, 7))
        ax = fig.add_subplot(111)
        numtins = 0
        patches = []
        lonlistarr = []
        latlistarr = []
        elevlistarr = []
        areaslistarr = []
        effareaslistarr = []
        xarr = []
        yarr = []
        zarr = []
        for c in self._pv:
            xarr.append(c[0])
            yarr.append(c[1])
            zarr.append(c[2])
        self._lon, self._lat = self.lonlat(xarr, yarr,projstring)
        lonlistarr.append(xarr)
        latlistarr.append(yarr)
        elevlistarr.append(zarr)
        area = []
        effarea = []
        zfarr = []
        for k in self._fv:
            v0 = self._pv[k[0]]
            v1 = self._pv[k[1]]
            v2 = self._pv[k[2]]
            self._arr.append(self.areaf(v0[0], v0[1], v0[2], v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]))
            self._effarr.append(self.areaf(v0[0], v0[1], 0.0, v1[0], v1[1], 0.0, v2[0], v2[1], 0.0))
            mid_point = np.asarray(
                    [(v0[0] + v1[0] + v2[0]) / 3, (v0[1] + v1[1] + v2[1]) / 3, (v0[2] + v1[2] + v2[2]) / 3])
            zfarr.append(mid_point[2])
        self._zfaces = np.asarray(zfarr)

        areaslistarr.append(np.asarray(area))
        effareaslistarr.append(np.asarray(effarea))
        #print(cv)
        #rgb = tuple(cv)
        rgb = self._materials
        #print(rgb)

        i = 0
        rgb_norm = []
        for (v, n, r, g, b) in land_covers:
            label = n
            rgb_1 = (r / 255, g / 255, b / 255)
            rgb_norm.append(rgb_1)
            patches.append(mpatches.Patch(color=rgb_1, label=label))
        i += 1

        cmap1 = colors.ListedColormap(rgb)
        cnames = []
        colors_floats = []
        for i in range(cmap1.N):
            single_rgb = (cmap1(i)[:3][0]/255, cmap1(i)[:3][1]/255, cmap1(i)[:3][2]/255)
            colors_floats.append(np.sqrt((cmap1(i)[:3][0]/255)**2+(cmap1(i)[:3][1]/255)**2+(cmap1(i)[:3][2]/255)**2))
            cnames.append(colors.rgb2hex(single_rgb))

        hex_list = cnames
        cmap = get_continuous_cmap(hex_list)

        plt.tripcolor(self._lon, self._lat, self._fv, facecolors=colors_floats, edgecolors='k', cmap=cmap)

        print(len(self._zfaces))
        print(len(self._fv))

        numtins += len(self._fv)
        title = 'Ncells: ' + str(numtins)
        plt.xlabel('Longitude, [deg]', fontsize=16)
        plt.ylabel('Latitude, [deg]', fontsize=16)
        plt.title(title, fontsize=18)
        #label = plabels[i]
        print(land_covers)
        ax.legend(handles=patches, loc=4)
        self._elevarr = np.array([elem for singleList in elevlistarr for elem in singleList])
        figname = jpg_folder + '/' + tid + '.jpg'
        #fig.savefig(figname)
        plt.show()

    def calculate_statistics(self):
        print("======== Effective area statistics ============")
        print(np.asarray(self._arr).mean())
        print(np.asarray(self._effarr).mean())
        print("mean: " + str(np.asarray(self._effarr).mean()))
        print("median: " + str(np.median(np.asarray(self._effarr))))
        print("variance: " + str(np.var(np.asarray(self._effarr))))
        print("standard deviation: " + str(np.std(np.asarray(self._effarr))))
        print("skew: " + str(stats.skew(np.asarray(self._effarr))))
        print("kurtosis: " + str(stats.kurtosis(np.asarray(self._effarr))))
        print("===== Normalised stats ======")
        neffarr = minmax_scale(np.asarray(self._effarr))
        print("mean: " + str(np.asarray(neffarr).mean()))
        print("median: " + str(np.median(np.asarray(neffarr))))
        print("variance: " + str(np.var(np.asarray(neffarr))))
        print("standard deviation: " + str(np.std(np.asarray(neffarr))))
        print("skew: " + str(stats.skew(np.asarray(neffarr))))
        print("kurtosis: " + str(stats.kurtosis(np.asarray(neffarr))))
        print("======== Effective height statistics ============")
        print("mean: " + str(np.mean(self._elevarr)))
        print("median: " + str(np.median(self._elevarr)))
        print("variance: " + str(np.var(self._elevarr)))
        print("standard deviation: " + str(np.std(self._elevarr)))
        print("skew: " + str(stats.skew(self._elevarr)))
        print("kurtosis: " + str(stats.kurtosis(self._elevarr)))
        print("===== Normalised stats ======")
        nelevarr = minmax_scale(np.asarray(self._elevarr))
        print("mean: " + str(np.asarray(nelevarr).mean()))
        print("median: " + str(np.median(np.asarray(nelevarr))))
        print("variance: " + str(np.var(np.asarray(nelevarr))))
        print("standard deviation: " + str(np.std(np.asarray(nelevarr))))
        print("skew: " + str(stats.skew(np.asarray(nelevarr))))
        print("kurtosis: " + str(stats.kurtosis(np.asarray(nelevarr))))

    ############################################################
    # This is an old peace. Ment to output shyft results on top of teh TIN figure
    # TODO: redo for teh new codebase bith shyft and tins
    # def view_tins_zfaces(self,tin_data_folder, tid, jpg_folder):
    #
    #     folder = '/Users/olgasilantyeva/[Documents]/Science/geohyd/projects/narayani/results/'
    #     figfolder2000 = '/Users/olgasilantyeva/[Documents]/Science/geohyd/projects/narayani/results/rptgsk/cid-10/figures-year-2000-spinned/'
    #     figfolder2003 = '/Users/olgasilantyeva/[Documents]/Science/geohyd/projects/narayani/results/rptgsk/cid-10/figures-year-2003-spinned/'
    #     figfolder = figfolder2000
    #
    #     file = 'rptgsk/cid-10/tin-slr-wfg-slope-asp180-2000-spinned/'
    #     csvfile = 'snow_spat_avg.csv'
    #     # csvfile = 'snow_spat_avg_m1.csv'
    #     csvfilerad = 'rad_spat_avg.csv'
    #     snow_spat_avg = pd.read_csv(folder + file + csvfile)
    #     rad_spat_avg = pd.read_csv(folder + file + csvfilerad)
    #     file1 = 'rptgsk/cid-10/tin-slr-wfg-aspect-2000-spinned/'
    #     snow_spat_avg1 = pd.read_csv(folder + file1 + csvfile)
    #
    #     tin_repo = tin.TinRepository(path=Path(tin_data_folder))
    #     geomdic = tin_repo.read(uid=tid)
    #     sns.set(font_scale=1.6)
    #     fig = plt.figure(figsize=(10, 7))
    #     ax = fig.add_subplot(111)
    #     i = 0
    #     numtins = 0
    #     plabels = []
    #     patches = []
    #     lonlistarr = []
    #     latlistarr = []
    #     elevlistarr = []
    #     areaslistarr = []
    #     effareaslistarr = []
    #     pointer = 0
    #     for k, v in geomdic.items():
    #         # print(k)
    #         plabels.append(str(k))
    #         pv = np.asarray(v.points)
    #         fv = np.asarray(v.faces)
    #         cv = np.asarray(v.colors[:][0])
    #         xarr = []
    #         yarr = []
    #         zarr = []
    #         for c in pv:
    #             xarr.append(c[0])
    #             yarr.append(c[1])
    #             zarr.append(c[2])
    #         lon, lat = lonlat(xarr, yarr)
    #         lonlistarr.append(xarr)
    #         latlistarr.append(yarr)
    #         elevlistarr.append(zarr)
    #         area = []
    #         effarea = []
    #         zfarr = []
    #         for k in fv:
    #             v0 = pv[k[0]]
    #             v1 = pv[k[1]]
    #             v2 = pv[k[2]]
    #             area.append(areaf(v0[0], v0[1], v0[2], v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]))
    #             effarea.append(areaf(v0[0], v0[1], 0.0, v1[0], v1[1], 0.0, v2[0], v2[1], 0.0))
    #             mid_point = np.asarray(
    #                 [(v0[0] + v1[0] + v2[0]) / 3, (v0[1] + v1[1] + v2[1]) / 3, (v0[2] + v1[2] + v2[2]) / 3])
    #             zfarr.append(mid_point[2])
    #         zfaces = np.asarray(zfarr)
    #         scaarr = []
    #         swearr = []
    #         swediffarr = []
    #         radtarr = []
    #         radparr = []
    #         for kk in range(len(zfaces)):
    #             scaarr.append(snow_spat_avg.snow_sca_spat_avg[pointer + kk])
    #             swearr.append(snow_spat_avg.snow_swe_spat_avg[pointer + kk])
    #             swediffarr.append(
    #                 abs(snow_spat_avg.snow_swe_spat_avg[pointer + kk] - snow_spat_avg1.snow_swe_spat_avg[pointer + kk]))
    #             radtarr.append(rad_spat_avg.rad_sim_t_spat_avg[pointer + kk])
    #             radparr.append(rad_spat_avg.rad_sim_cs_p_spat_avg[pointer + kk])
    #         pointer += len(zfaces)
    #         sca = np.asarray(scaarr)
    #         swe = np.asarray(swearr)
    #         radt = np.asarray(radtarr)
    #         radp = np.asarray(radparr)
    #         swediff = np.asarray(swediffarr)
    #         print(len(swe))
    #         print(len(zfaces))
    #         print(swe)
    #         # print(zfaces)
    #         areaslistarr.append(np.asarray(area))
    #         effareaslistarr.append(np.asarray(effarea))
    #         rgb = tuple(cv)
    #         norm = colors.Normalize(vmin=0, vmax=7000)
    #         cmap = colors.ListedColormap([rgb])
    #         fc = np.full(len(fv), 1)
    #         print(len(zfaces))
    #         print('max swe: ' + str(max(swe)))
    #         print('min swe: ' + str(min(swe)))
    #         print('max swediff: ' + str(max(swediff)))
    #         print('min swediff: ' + str(min(swediff)))
    #         # plt.tripcolor(lon, lat, fv, facecolors=swe, edgecolors='k', cmap='tab20b',vmin=0, vmax = 80)
    #         # sca
    #         # plt.tripcolor(lon, lat, fv, facecolors=sca, edgecolors='k', cmap='tab20b', vmin=0, vmax=1)
    #         # plt.tripcolor(lon, lat, fv, facecolors=swediff, edgecolors='k', cmap='tab20b', vmin=0, vmax=6)
    #         # plt.tripcolor(lon, lat, fv, facecolors=radt, edgecolors='k', cmap='coolwarm', vmin=100, vmax=260)
    #         plt.tripcolor(lon, lat, fv, facecolors=radp, edgecolors='k', cmap='coolwarm', vmin=100, vmax=350)
    #         numtins += len(fv)
    #         title = 'Ncells: ' + str(numtins)
    #         plt.xlabel('Longitude, [deg]', fontsize=16)
    #         plt.ylabel('Latitude, [deg]', fontsize=16)
    #         plt.title(title, fontsize=18)
    #         # label = plabels[i]
    #         # patches.append(mpatches.Patch(color=rgb, label=label))
    #         i += 1
    #
    #     # ax.legend(handles=patches, loc=3)
    #     plt.colorbar()
    #     # figname = figfolder + 'tin-mesh-snow-sca-slr.jpg'
    #     # figname = figfolder + 'tin-mesh-rad_t-slr-slope-asp180.jpg'
    #     figname = figfolder + 'tin-mesh-rad_p-slr-slope-asp180.jpg'
    #     fig.savefig(figname)