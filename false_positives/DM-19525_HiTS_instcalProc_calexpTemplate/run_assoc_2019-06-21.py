import lsst.daf.persistence as dafPersist
import lsst.daf.persistence.butlerExceptions as butlerExceptions
from lsst.ap.association import (
    AssociationTask,
    MapDiaSourceTask,
    make_dia_object_schema,
    make_dia_source_schema)
from lsst.dax.ppdb import Ppdb, PpdbConfig
from lsst.utils import getPackageDir
from lsst.pipe.tasks.imageDifference import (ImageDifferenceTask,
                                             ImageDifferenceConfig)
import lsst.log as Log
import os


def _data_file_name(basename, module_name):
    """Return path name of a data file.

    Parameters
    ----------
    basename : `str`
        Name of the file to add to the path string.
    module_name : `str`
        Name of lsst stack package environment variable.

    Returns
    -------
    data_file_path : `str`
       Full path of the file to load from the "data" directory in a given
       repository.
    """
    return os.path.join(getPackageDir(module_name), "data", basename)


class RunAssociation:

    def __init__(self, db_file):
        """Create similar configuration for tasks ad in ap_pipe.
        """

        self.log = Log.getLogger("RunAssociation")
        self.ppdbConfig = PpdbConfig()
        self.ppdbConfig.db_url = "sqlite:///" + db_file
        self.ppdbConfig.isolation_level = "READ_UNCOMMITTED"
        self.ppdbConfig.dia_object_index = "baseline"
        self.ppdbConfig.dia_object_columns = []
        self.ppdbConfig.connection_timeout = 240
        self.ppdbConfig.schema_file = _data_file_name(
            "ppdb-schema.yaml", "dax_ppdb")
        self.ppdbConfig.column_map = _data_file_name(
            "ppdb-ap-pipe-afw-map.yaml", "ap_association")
        self.ppdbConfig.extra_schema_file = _data_file_name(
            "ppdb-ap-pipe-schema-extra.yaml", "ap_association")

        self.ppdb = Ppdb(config=self.ppdbConfig,
                         afw_schemas=dict(DiaObject=make_dia_object_schema(),
                                          DiaSource=make_dia_source_schema()))
        # ppdb.makeSchema()
        self.differencerConfig = ImageDifferenceConfig()
        # Schema is different if we do decorrelation
        self.differencerConfig.doDecorrelation = True
        self.differencerSchema = ImageDifferenceTask(config=self.differencerConfig).schema
        self.diaSourceDpddifier = MapDiaSourceTask(inputSchema=self.differencerSchema)
        self.associator = AssociationTask()

        self.diffType = "deep"

    def procOneExposure(self, dataId):
        """Run AssociationTask on the exposure selected by ``dataId``.
        """

        try:
            catalog = self.butler.get(self.diffType + "Diff_diaSrc", dataId=dataId)
            diffim = self.butler.get(self.diffType + "Diff_differenceExp", dataId=dataId)
        except butlerExceptions.NoResults:
            self.log.info("Data does not exist %s", dataId)
            return

        dia_sources = self.diaSourceDpddifier.run(catalog,
                                                  diffim,
                                                  return_pandas=True)
        self.associator.run(dia_sources, diffim, self.ppdb)

    def run(self, visitIds):
        """Loop through visits and process all available exposures within.
        """
        self.butler = dafPersist.Butler("repo/rerun/imgDiff_2019-06-20")
        for v in visitIds:
            ids = self.butler.queryMetadata(self.diffType + "Diff_diaSrc",
                                            ['visit', 'ccdnum', 'filter'],
                                            {'visit': v})
            for x in ids:
                dataId = {'visit': x[0], 'ccdnum': x[1], 'filter': x[2]}
                self.log.info("Processing %s", dataId)
                self.procOneExposure(dataId)


def main():
    """Run AssociationTask for the repository.

    Run from /project/gkovacs/hits_instcals.

    """
    # Log.getLogger('').setLevel(Log.DEBUG)
    Log.getLogger('association').setLevel(Log.DEBUG)
    assoc = RunAssociation("repo/rerun/imgDiff_2019-06-20/association.db")
    assoc.run((411420, 419802, 411371))


if __name__ == "__main__":
    main()
