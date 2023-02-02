# This file is part of Shyft. Copyright 2015-2018 SiH, JFB, OS, YAS, Statkraft AS
# See file COPYING for more details **/
"""
Read region rasputin files with cell data.

"""

from os import path
import numpy as np
from netCDF4 import Dataset
from shyft.hydrology.repository import interfaces
from pathlib import Path
from typing import Union
import pyproj
from functools import partial
from shapely.geometry import Polygon, MultiPolygon, box
from shapely.ops import transform

from shyft.hydrology import shyftdata_dir
from shyft.hydrology.orchestration.configuration.config_interfaces import RegionConfig, ModelConfig, RegionConfigError
from shyft.hydrology.orchestration.configuration.dict_configs import DictModelConfig, DictRegionConfig
from shyft.hydrology.repository.netcdf.utils import create_ncfile, make_proj


def get_land_type(rgb, land_covers_available)->(int,str):
    for (v,n,r,g,b) in land_covers_available:
        rgbs = (r,g, b)
        (value,name) = (v,n)
        if rgb == rgbs:
            return (value,name)
        else:
            pass


class CFRegionModelRepositoryError(Exception):
    pass


class CFRegionModelRepository(interfaces.RegionModelRepository):
    """
    Repository that delivers fully specified shyft api region_models
    based on data found in netcdf files.
    """

    def __init__(self, region, model):
        """
        Parameters
        ----------
        region: either a dictionary suitable to be instantiated as a
            RegionConfig Object or a sublcass of the interface RegionConfig
            containing regional information, like
            catchment overrides, and which netcdf file to read
        model: either a dictionary suitable to be instantiated as a
            ModelConfig Object or a subclass of interface ModelConfig
            Object containing model information, i.e.
            information concerning interpolation and model
            parameters
        """

        if not isinstance(region, RegionConfig):
            region_config = DictRegionConfig(region)
        if not isinstance(model, ModelConfig):
            model_config = DictModelConfig(model)
        else:
            region_config = region
            model_config = model

        if not isinstance(region_config, RegionConfig) or \
                not isinstance(model_config, ModelConfig):
            raise interfaces.InterfaceError()
        self._rconf = region_config
        self._mconf = model_config
        self._region_model = model_config.model_type()  # region_model
        self._mask = None
        self._epsg = self._rconf.domain()["EPSG"]  # epsg
        self._tin_uid_ = self._rconf.repository()["params"]["tin_uid"]
        filename = self._rconf.repository()["params"]["data_file"]
        filename = path.expandvars(filename)
        if not path.isabs(filename):
            # Relative paths will be prepended the data_dir
            filename = path.join(shyftdata_dir, filename)
        if not path.isfile(filename):
            raise CFRegionModelRepositoryError("No such file '{}'".format(filename))
        tinfolder = self._rconf.repository()["params"]["tin_data_folder"]
        tinfolder = path.expandvars(tinfolder)
        if not path.isabs(tinfolder):
            # Relative paths will be prepended the data_dir
            tinfolder = path.join(shyftdata_dir, tinfolder)
        if not path.isdir(tinfolder):
            raise CFRegionModelRepositoryError("No such file folder '{}'".format(tinfolder))
        self._data_file = filename
        self._tin_data_folder = tinfolder
        self._catch_ids = self._rconf.catchments()
        self._get_from_tin_ = self._rconf.repository()["params"]["get_model_from_tin_repo"]
        self.bounding_box = None

    def _limit(self, x, y, data_cs, target_cs):
        """
        Parameters
        ----------
        """
        # Get coordinate system for arome data
        data_proj = make_proj(data_cs)
        target_proj = make_proj(target_cs)

        # Find bounding box in arome projection
        bbox = self.bounding_box
        bb_proj = pyproj.transform(target_proj, data_proj, bbox[0], bbox[1])
        x_min, x_max = min(bb_proj[0]), max(bb_proj[0])
        y_min, y_max = min(bb_proj[1]), max(bb_proj[1])

        # Limit data
        xy_mask = ((x <= x_max) & (x >= x_min) & (y <= y_max) & (y >= y_min))

        xy_inds = np.nonzero(xy_mask)[0]

        # Transform from source coordinates to target coordinates
        xx, yy = pyproj.transform(data_proj, target_proj, x, y)

        return xx, yy, xy_mask, xy_inds

    def get_region_model(self, region_id, catchments=None):
        """
        Return a fully specified shyft api region_model for region_id, based on data found
        in netcdf dataset.

        Parameters
        -----------
        region_id: string
            unique identifier of region in data

        catchments: list of unique integers
            catchment indices when extracting a region consisting of a subset
            of the catchments has attribs to construct params and cells etc.

        Returns
        -------
        region_model: shyft.api type
        """

        if not (self._get_from_tin_):
            with Dataset(self._data_file) as dset:
                Vars = dset.variables
                c_ids = Vars["catchment_id"][:]
                xcoord = Vars['x'][:]
                ycoord = Vars['y'][:]
                m_catch = np.ones(len(c_ids), dtype=bool)
                if self._catch_ids is not None:
                    m_catch = np.in1d(c_ids, self._catch_ids)
                    xcoord_m = xcoord[m_catch]
                    ycoord_m = ycoord[m_catch]

                dataset_epsg = None
                if 'crs' in Vars.keys():
                    dataset_epsg = Vars['crs'].epsg_code.split(':')[1]
                if not dataset_epsg:
                    raise interfaces.InterfaceError("netcdf: epsg attr not found in group elevation")

                target_cs = f"EPSG:{self._epsg}"
                source_cs = f"EPSG:{dataset_epsg}"

                # Construct bounding region
                box_fields = set(("lower_left_x", "lower_left_y", "step_x", "step_y", "nx", "ny", "EPSG"))
                if box_fields.issubset(self._rconf.domain()):
                    tmp = self._rconf.domain()
                    epsg = tmp["EPSG"]
                    x_min = tmp["lower_left_x"]
                    x_max = x_min + tmp["nx"]*tmp["step_x"]
                    y_min = tmp["lower_left_y"]
                    y_max = y_min + tmp["ny"]*tmp["step_y"]
                    bounding_region = BoundingBoxRegion(np.array([x_min, x_max]),
                                                        np.array([y_min, y_max]), epsg, self._epsg)
                else:
                    bounding_region = BoundingBoxRegion(xcoord_m, ycoord_m, dataset_epsg, self._epsg)
                self.bounding_box = bounding_region.bounding_box(self._epsg)
                x, y, m_xy, _ = self._limit(xcoord, ycoord, source_cs, target_cs)
                mask = ((m_xy) & (m_catch))

                tin_x0_arr = Vars['x0'][mask]
                tin_y0_arr = Vars['y0'][mask]
                tin_z0_arr = Vars['z0'][mask]
                tin_x1_arr = Vars['x1'][mask]
                tin_y1_arr = Vars['y1'][mask]
                tin_z1_arr = Vars['z1'][mask]
                tin_x2_arr = Vars['x2'][mask]
                tin_y2_arr = Vars['y2'][mask]
                tin_z2_arr = Vars['z2'][mask]

                c_ids = Vars["catchment_id"][mask]
                c_ids_unique = list(np.unique(c_ids))
                # c_indx = np.array([c_ids_unique.index(cid) for cid in c_ids]) # since ID to Index conversion not necessary

                ff = Vars["forest-fraction"][mask]
                lf = Vars["lake-fraction"][mask]
                rf = Vars["reservoir-fraction"][mask]
                gf = Vars["glacier-fraction"][mask]
                uf = np.full(len(ff), 0)
        else:
            tin_x0_arr, tin_y0_arr, tin_z0_arr, tin_x1_arr, tin_y1_arr, tin_z1_arr, tin_x2_arr, tin_y2_arr, tin_z2_arr, lf, gf, rf, ff, uf, c_ids, c_ids_unique, bounding_region = self.parse_tinrepo()

        # Construct region parameter:
        region_parameter = self._region_model.parameter_t()
        for p_type_name, value_ in self._mconf.model_parameters().items():
            if hasattr(region_parameter, p_type_name):
                sub_param = getattr(region_parameter, p_type_name)
                for p, v in value_.items():
                    if hasattr(sub_param, p):
                        setattr(sub_param, p, v)
                    else:
                        raise RegionConfigError("Invalid parameter '{}' for parameter set '{}'".format(p, p_type_name))
            else:
                raise RegionConfigError("Invalid parameter set '{}' for selected model '{}'".format(p_type_name, self._region_model.__name__))

        radiation_slope_factor = 1.0

        cell_geo_data = np.column_stack([tin_x0_arr, tin_y0_arr, tin_z0_arr, tin_x1_arr, tin_y1_arr, tin_z1_arr, tin_x2_arr, tin_y2_arr, tin_z2_arr, np.full(len(tin_x0_arr), self._epsg), c_ids, gf, lf, rf, ff])

        cell_vector = self._region_model.cell_t.vector_t.create_from_geo_cell_data_vector_to_tin(np.ravel(cell_geo_data))

        # Construct catchment overrides
        catchment_parameters = self._region_model.parameter_t.map_t()
        for cid, catch_param in self._rconf.parameter_overrides().items():
            if cid in c_ids_unique:
                param = self._region_model.parameter_t(region_parameter)
                for p_type_name, value_ in catch_param.items():
                    if hasattr(param, p_type_name):
                        sub_param = getattr(param, p_type_name)
                        for p, v in value_.items():
                            if hasattr(sub_param, p):
                                setattr(sub_param, p, v)
                            else:
                                raise RegionConfigError("Invalid parameter '{}' for catchment parameter set '{}'".format(p, p_type_name))
                    else:
                        raise RegionConfigError("Invalid catchment parameter set '{}' for selected model '{}'".format(p_type_name, self._region_model.__name__))

                catchment_parameters[cid] = param
        region_model = self._region_model(cell_vector, region_parameter, catchment_parameters)
        region_model.bounding_region = bounding_region
        region_model.catchment_id_map = c_ids_unique

        def do_clone(x):
            clone = x.__class__(x)
            clone.bounding_region = x.bounding_region
            clone.catchment_id_map = x.catchment_id_map
            # clone.gis_info = polygons  # cell shapes not included yet
            return clone

        region_model.clone = do_clone
        return region_model

    def cell_data_to_netcdf(self, region_model, output_dir):
        """
        Writes cell_data from a shyft region_model in the same format the
         'cf_region_model_repository' expects.

        Parameters
        -----------
        region_model: shyft.region_model

        model_id: str identifier of region_model

        Returns
        -------


        """

        nc_file = "{}_cell_data.nc".format(output_dir)
        # repository = {'class': self.__class__,
        #               'params': {'data_file': nc_file}
        #               }

        # with open('{}.yaml'.format(nc_file), 'w') as yml:
        #     yaml.dump(repository, yml)

        dimensions = {'cell': len(region_model.cells)}

        variables = {'y': [np.float, ('cell',), {'axis': 'mid_Y',
                                                 'units': 'm',
                                                 'standard_name': 'projection_y_coordinate'}],

                     'x': [np.float, ('cell',), {'axis': 'mid_X',
                                                 'units': 'm',
                                                 'standard_name': 'projection_x_coordinate'}],

                     'z': [np.float, ('cell',), {'axis': 'mid_Z',
                                                 'units': 'm',
                                                 'standard_name': 'height',
                                                 'long_name': 'height above mean sea level'}],

                     'crs': [np.int32, ('cell',), {'grid_mapping_name': 'transverse_mercator',
                                                   'epsg_code': 'EPSG:' + str(region_model.bounding_region.epsg()),
                                                   'proj4': "+proj = utm + zone = 33 + datum = WGS84 + units = m + no_defs + ellps = WGS84 + towgs84=0,0,0"}],

                     'area': [np.float, ('cell',), {'grid_mapping': 'crs',
                                                    'units': 'm^2',
                                                    'coordinates': 'y x z'}],

                     'forest-fraction': [np.float, ('cell',), {'grid_mapping': 'crs',
                                                               'units': '-',
                                                               'coordinates': 'y x z'}],

                     'glacier-fraction': [np.float, ('cell',), {'grid_mapping': 'crs',
                                                                'units': '-',
                                                                'coordinates': 'y x z'}],

                     'lake-fraction': [np.float, ('cell',), {'grid_mapping': 'crs',
                                                             'units': '-',
                                                             'coordinates': 'y x z'}],

                     'reservoir-fraction': [np.float, ('cell',), {'grid_mapping': 'crs',
                                                                  'units': '-',
                                                                  'coordinates': 'y x z'}],

                     'catchment_id': [np.int32, ('cell',), {'grid_mapping': 'crs',
                                                            'units': '-',
                                                            'coordinates': 'y x z'}],
                     'x0': [np.float, ('cell',), {'axis': 'x0',
                                                  'units': 'm',
                                                  'standard_name': 'x0'}],

                     'y0': [np.float, ('cell',), {'axis': 'y0',
                                                  'units': 'm',
                                                  'standard_name': 'y0'}],

                     'z0': [np.float, ('cell',), {'axis': 'z0',
                                                  'units': 'm',
                                                  'standard_name': 'height0',
                                                  'long_name': 'z0'}],
                     'x1': [np.float, ('cell',), {'axis': 'x1',
                                                  'units': 'm',
                                                  'standard_name': 'x1'}],

                     'y1': [np.float, ('cell',), {'axis': 'y1',
                                                  'units': 'm',
                                                  'standard_name': 'y1'}],

                     'z1': [np.float, ('cell',), {'axis': 'z1',
                                                  'units': 'm',
                                                  'standard_name': 'height1',
                                                  'long_name': 'z1'}],
                     'x2': [np.float, ('cell',), {'axis': 'x2',
                                                  'units': 'm',
                                                  'standard_name': 'x2'}],

                     'y2': [np.float, ('cell',), {'axis': 'y2',
                                                  'units': 'm',
                                                  'standard_name': 'y2'}],

                     'z2': [np.float, ('cell',), {'axis': 'z2',
                                                  'units': 'm',
                                                  'standard_name': 'height2',
                                                  'long_name': 'z2'}],
                     'slopes': [np.float, ('cell',), {'units': 'deg',
                                                      'standard_name': 'slopes'}],
                     'aspects': [np.float, ('cell',), {'units': 'deg',
                                                       'standard_name': 'aspects'}],
                     }

        create_ncfile(nc_file, variables, dimensions)
        nci = Dataset(nc_file, 'a')

        extracted_geo_cell_data = region_model.extract_geo_cell_data()
        nci.variables['x'][:] = [gcd.mid_point().x for gcd in extracted_geo_cell_data]
        nci.variables['y'][:] = [gcd.mid_point().y for gcd in extracted_geo_cell_data]
        nci.variables['z'][:] = [gcd.mid_point().z for gcd in extracted_geo_cell_data]
        nci.variables['area'][:] = [gcd.area() for gcd in extracted_geo_cell_data]
        nci.variables['catchment_id'][:] = [gcd.catchment_id() for gcd in extracted_geo_cell_data]
        nci.variables['lake-fraction'][:] = [gcd.land_type_fractions_info().lake() for gcd in extracted_geo_cell_data]
        nci.variables['reservoir-fraction'][:] = [gcd.land_type_fractions_info().reservoir() for gcd in extracted_geo_cell_data]
        nci.variables['glacier-fraction'][:] = [gcd.land_type_fractions_info().glacier() for gcd in extracted_geo_cell_data]
        nci.variables['forest-fraction'][:] = [gcd.land_type_fractions_info().forest() for gcd in extracted_geo_cell_data]
        nci.variables['x0'][:] = [gcd.vertexes()[0].x for gcd in extracted_geo_cell_data]
        nci.variables['y0'][:] = [gcd.vertexes()[0].y for gcd in extracted_geo_cell_data]
        nci.variables['z0'][:] = [gcd.vertexes()[0].z for gcd in extracted_geo_cell_data]
        nci.variables['x1'][:] = [gcd.vertexes()[1].x for gcd in extracted_geo_cell_data]
        nci.variables['y1'][:] = [gcd.vertexes()[1].y for gcd in extracted_geo_cell_data]
        nci.variables['z1'][:] = [gcd.vertexes()[1].z for gcd in extracted_geo_cell_data]
        nci.variables['x2'][:] = [gcd.vertexes()[2].x for gcd in extracted_geo_cell_data]
        nci.variables['y2'][:] = [gcd.vertexes()[2].y for gcd in extracted_geo_cell_data]
        nci.variables['z2'][:] = [gcd.vertexes()[2].z for gcd in extracted_geo_cell_data]
        nci.variables['slopes'][:] = [gcd.slope() for gcd in extracted_geo_cell_data]
        nci.variables['aspects'][:] = [gcd.aspect() for gcd in extracted_geo_cell_data]

        nci.close()

    def parse_tinrepo(self):
        # instead of rasputin we use pip available h5py to parse h5 files
        import h5py
        # Construct cells from TIN repository
        # class GlobCovLandTypes(Enum):
        #     crop_type_1 = 11
        #     crop_type_2 = 14
        #     crop_type_3 = 20
        #     crop_type_4 = 30
        #     forest_type_1 = 40
        #     forest_type_2 = 50
        #     forest_type_3 = 60
        #     forest_type_4 = 70
        #     forest_type_5 = 90
        #     forest_type_6 = 100
        #     shrub_type_1 = 110
        #     shrub_type_2 = 120
        #     shrub_type_3 = 130
        #     vegetation_type_1 = 140
        #     vegetation_type_2 = 150
        #     flood_type_1 = 160
        #     flood_type_2 = 170
        #     flood_type_3 = 180
        #     artificial = 190
        #     bare = 200
        #     water = 210
        #     snow_and_ice = 220
        #     no_data = 230

        # class CorineLandCoverType():
        #     # Artificial
        #     urban_fabric_cont = 111
        #     urban_fabric_discont = 112
        #     industrial_unit = 121
        #     road_and_rail = 122
        #     port = 123
        #     airport = 124
        #     mineral_extraction = 131
        #     dump_site = 132
        #     constrution_site = 133
        #     urban_green = 141
        #     sport_and_leisure = 142
        #
        #     # Agricultural
        #     arable_land_non_irr = 211
        #     permanent_irr = 212
        #     rice_field = 213
        #     vinyard = 221
        #     fruit_and_berry = 222
        #     olive_grove = 223
        #     pasture = 231
        #     mix_annual_permament_crop = 241
        #     complex_cultivation = 242
        #     mix_agri_natural = 243
        #     agro_forestry = 244
        #
        #     # Forest
        #     broad_leaved = 311
        #     coniferous = 312
        #     mixed_forest = 313
        #     natural_grass = 321
        #     moors_and_heath = 322
        #     sclerophyllous = 323
        #     transitional_woodland_shrub = 324
        #     beach_dune_sand = 331
        #     bare_rock = 332
        #     sparse_veg = 333
        #     burnt = 334
        #     glacier_and_snow = 335
        #
        #     # Wetland
        #     inland_march = 411
        #     peat_bog = 412
        #     salt_march = 421
        #     saline = 422
        #     intertidal_flat = 423
        #
        #     # Water bodies
        #     water_course = 511
        #     water_body = 512
        #     coastal_lagoon = 521
        #     estuary = 522
        #     sea_and_ocean = 523

        # TODO: define more types on c++ side and make proper conversion from orchestration
        # Right now Rasputin is not able to delineate watersheds, so the solution is to put all subcatchments insame directory and read one-by-one
        # cid = tin_uid = filename --> key to find subcatchment

        tin_v0_arr = []
        tin_v1_arr = []
        tin_v2_arr = []
        lf = []
        rf = []
        ff = []
        gf = []
        uf = []
        c_ids = []
        #print(self._catch_ids)
        for tid in self._tin_uid_:
            tin_repo = h5py.File(Path(self._tin_data_folder +"/"+ tid + ".h5"), 'r')
            # from here I already know teh structure of the file, by asking tin_repo.keys() one can check the main groups
            # in the rasputin generated tins main groups are "information" and "tin"
            tins = tin_repo.get("tin")
            # tins.keys() return "face_fields", "faces", "points"
            points = tins.get("points")  # this is already a dataset, one cah check points.shape,
            # in the rasputin generated tin the shape is normally (L,3), where L -- is number of vertexes
            xarr = points[:, 0]
            yarr = points[:, 1]
            zarr = points[:, 2]
            faces = tins.get('faces')
            # now we got all our vertexes, but still need to gain information about land_types
            # colors and codes for land_cover type are kept with tins
            face_fields = tins.get(
                "face_fields")  # this one carries additional information about cover_color and cover_type
            tincolors = face_fields.get("cover_color")  # shape is same as for faces, as it is a color of face
            land_cover_type = face_fields.get("cover_type")
            # but text information about the name of the land_type is kept in information group
            info = tin_repo.get("information")
            # info.keys() returns "land_covers"
            land_covers = info.get("land_covers")

            materials = [(int(np.round(r * 255)), int(np.round(g * 255)), int(np.round(b * 255))) for (r, g, b) in
                         tincolors]
            cid_val = int(''.join(filter(str.isdigit, tid)))
            pv = np.asarray(points)
            fv = np.asarray(faces)
            #projstring = "EPSG:32645"  # default to Nepal
            projstring = "zone=45"
            land_covers_names = []
            for c in materials:
                (value, name) = get_land_type(c, land_covers)
                land_covers_names.append(name.decode('utf-8'))
                print(name.decode('utf-8'))

            # source_cs = v.projection
            # TODO: if rasputin ends u with epsg string change here:
            # source_cs = "EPSG:32633"
            # crs = pyproj.Proj(source_cs, preserve_units=True)
            #dataset_epsg = projstring
            dataset_epsg = None
            # print(projstring)
            spstr = projstring.split(" ")
            for k in spstr:
                if ('zone=' in k):
                    zone = int(k.split('=')[1])
            # print(crs)
            dataset_epsg = 32600 + zone
            # print(dataset_epsg)
            if not dataset_epsg:
                raise interfaces.InterfaceError("netcdf: can't define epsg from proj string")

            pp = np.asarray(pv)[np.asarray(fv).reshape(-1)].reshape(-1, 3)
            v0 = np.asarray(pp[0::3]).reshape(-1).reshape(-1, 3)
            v1 = np.asarray(pp[1::3]).reshape(-1).reshape(-1, 3)
            # print(v1)
            v2 = np.asarray(pp[2::3]).reshape(-1).reshape(-1, 3)
            tin_v0_arr = np.append(tin_v0_arr, v0)
            tin_v1_arr = np.append(tin_v1_arr, v1)
            tin_v2_arr = np.append(tin_v2_arr, v2)
            # TODO: make more land cover types on shyft side in main branch
            # parsing to shyft internal types: lf -lake/water body, ff - forest, rf - reservoir, gf - glacier, uf - undefiened
            for lc in land_covers_names:
                print(lc)
                if lc=="snow_and_ice":
                    lf = np.append(lf, 0)
                    ff = np.append(ff, 0)
                    rf = np.append(rf, 0)
                    gf = np.append(gf, 1)
                    uf = np.append(uf, 0)
                elif lc=="bare":
                    lf = np.append(lf, 0)
                    ff = np.append(ff, 0)
                    rf = np.append(rf, 0)
                    gf = np.append(gf, 0)
                    uf = np.append(uf, 1)
                else:
                    lf = np.append(lf, 0)
                    ff = np.append(ff, 1)
                    rf = np.append(rf, 0)
                    gf = np.append(gf, 0)
                    uf = np.append(uf, 0)
            print(ff)
            c_ids = np.append(c_ids,np.full(len(v0),cid_val))

            tin_x0_arr = tin_v0_arr.reshape(-1)[0::3]
            tin_y0_arr = tin_v0_arr.reshape(-1)[1::3]
            tin_z0_arr = tin_v0_arr.reshape(-1)[2::3]
            tin_x1_arr = tin_v1_arr.reshape(-1)[0::3]
            tin_y1_arr = tin_v1_arr.reshape(-1)[1::3]
            tin_z1_arr = tin_v1_arr.reshape(-1)[2::3]
            tin_x2_arr = tin_v2_arr.reshape(-1)[0::3]
            tin_y2_arr = tin_v2_arr.reshape(-1)[1::3]
            tin_z2_arr = tin_v2_arr.reshape(-1)[2::3]
            c_ids_unique = list(np.unique(c_ids))
            target_cs = "EPSG:{}".format(self._epsg)
            source_cs = "EPSG:{}".format(dataset_epsg)
            self._epsg = dataset_epsg    
        # Construct bounding region
        c_ids = c_ids.reshape(-1)
        m_catch = np.ones(len(c_ids), dtype=bool)
        box_fields = set(("lower_left_x", "lower_left_y", "step_x", "step_y", "nx", "ny", "EPSG"))
        if box_fields.issubset(self._rconf.domain()):
            tmp = self._rconf.domain()
            epsg = tmp["EPSG"]
            x_min = tmp["lower_left_x"]
            x_max = x_min + tmp["nx"]*tmp["step_x"]
            y_min = tmp["lower_left_y"]
            y_max = y_min + tmp["ny"]*tmp["step_y"]
            bounding_region = BoundingBoxRegion(np.array([x_min, x_max]),
                                                np.array([y_min, y_max]), epsg, self._epsg)
        else:
            bounding_region = BoundingBoxRegion(tin_x0_arr, tin_y0_arr, dataset_epsg, self._epsg)
        self.bounding_box = bounding_region.bounding_box(self._epsg)
        x, y, m_xy, _ = self._limit(tin_x0_arr, tin_y0_arr, source_cs, target_cs)
        mask = ((m_xy) & (m_catch))  # TODO utilize mask to run only on cells included into bb
        return tin_x0_arr, tin_y0_arr, tin_z0_arr, tin_x1_arr, tin_y1_arr, tin_z1_arr, tin_x2_arr, tin_y2_arr, tin_z2_arr, lf, gf, rf, ff, uf, c_ids, c_ids_unique, bounding_region


class BoundingBoxRegion(interfaces.BoundingRegion):

    def __init__(self, x, y, point_epsg, target_epsg):
        self._epsg = str(point_epsg)
        x_min = x.ravel().min()
        x_max = x.ravel().max()
        y_min = y.ravel().min()
        y_max = y.ravel().max()
        self.x = np.array([x_min, x_max, x_max, x_min], dtype="d")
        self.y = np.array([y_min, y_min, y_max, y_max], dtype="d")
        self.x, self.y = self.bounding_box(target_epsg)
        self._epsg = str(target_epsg)
        self._polygon = box(x_min, y_min, x_max, y_max)

    def bounding_box(self, epsg):
        epsg = str(epsg)
        if epsg == self.epsg():
            return np.array(self.x), np.array(self.y)
        else:
            source_proj = make_proj(f"EPSG:{self.epsg()}")
            target_proj = make_proj(f"EPSG:{epsg}")
            return [np.array(a) for a in pyproj.transform(source_proj, target_proj, self.x, self.y)]

    def bounding_polygon(self, epsg: int) -> Union[Polygon, MultiPolygon]:
        """Implementation of interface.BoundingRegion"""
        if epsg == self.epsg():
            return self._polygon
        else:
            source_proj = make_proj(f"EPSG:{self.epsg()}")
            target_proj = make_proj(f"EPSG:{epsg}")
            project = partial(pyproj.transform, source_proj, target_proj)
            return transform(project, self._polygon)

    def epsg(self):
        return self._epsg
