from peewee import *

database = MySQLDatabase('mshakya', **{'passwd': 'f001r60', 'host': 'localhost', 'user': 'mshakya'})

class UnknownField(object):
    pass

class BaseModel(Model):
    class Meta:
        database = database

class Accgitable(BaseModel):
    accid = CharField(db_column='accID', max_length=15, null=True)
    gi = IntegerField(null=True)

    class Meta:
        db_table = 'AccGiTable'

class Accpostable(BaseModel):
    accid = CharField(db_column='accID', max_length=15, null=True)
    position = IntegerField(null=True)

    class Meta:
        db_table = 'AccPosTable'

class Acctaxidtable(BaseModel):
    accid = CharField(db_column='accID', max_length=15, null=True)
    taxid = IntegerField(db_column='taxID', null=True)

    class Meta:
        db_table = 'AccTaxIDTable'

class Acctitletable(BaseModel):
    title = TextField(db_column='Title', null=True)
    accid = CharField(db_column='accID', max_length=15, null=True)

    class Meta:
        db_table = 'AccTitleTable'

class Accessiongitable(BaseModel):
    accessionid = TextField(db_column='accessionID', null=True)
    gi = IntegerField(null=True)

    class Meta:
        db_table = 'AccessionGiTable'

class Clustertable(BaseModel):
    position = IntegerField(db_column='Position', null=True)
    accid = CharField(db_column='accID', max_length=15, null=True)
    cluster = IntegerField(null=True)
    gi = PrimaryKeyField()
    homolog = TextField(null=True)
    method = TextField(null=True)
    org2 = TextField(null=True)
    organism = TextField(null=True)
    pos2 = IntegerField(null=True)
    status = TextField(null=True)

    class Meta:
        db_table = 'ClusterTable'

class Gtasfinalmanual(BaseModel):
    family = TextField(db_column='Family', null=True)
    genus = TextField(db_column='Genus', null=True)
    orders = TextField(db_column='Orders', null=True)
    end = IntegerField(null=True)
    gi = IntegerField(null=True)
    homolog = TextField(null=True)
    method = TextField(null=True)
    organism = TextField(null=True)
    start = IntegerField(null=True)
    status = TextField(null=True)
    strands = CharField(max_length=1, null=True)

    class Meta:
        db_table = 'GTAsFinalManual'

class Gihomolog(BaseModel):
    gi = IntegerField(db_column='gi')
    homolog = TextField(null=True)
    method = TextField(null=True)
    organism = TextField(null=True)
    status = TextField(null=True)

    class Meta:
        db_table = 'GiHomolog'

class Giorganism(BaseModel):
    gi = PrimaryKeyField()
    organism = TextField()
    status = TextField()

    class Meta:
        db_table = 'GiOrganism'

class Mlpredictions(BaseModel):
    thous_best_bns_selected_feats = CharField(db_column='1000_best_BNS_selected_feats', max_length=10, null=True)
    confidence_1000_feats = FloatField(db_column='Confidence_1000_feats', null=True)
    all_feats_minus_sparse = CharField(max_length=10, null=True)
    gi = IntegerField(primary_key=True, null=True, db_column='gid')
    homolog = CharField(max_length=7, null=True)
    cluster_size = IntegerField(db_column='cluster_size')

    class Meta:
        db_table = 'MLPredictions'

class Ptttables(BaseModel):
    cog = TextField(db_column='COG', null=True)
    function = TextField(db_column='Function', null=True)
    position = IntegerField(db_column='Position', null=True)
    product = TextField(db_column='Product', null=True)
    synonymcode = TextField(db_column='SynonymCode', null=True)
    end = IntegerField(null=True)
    gene = TextField(null=True)
    gi = PrimaryKeyField(db_column='gi')
    length = IntegerField(null=True)
    start = IntegerField(null=True)
    strand = CharField(max_length=1, null=True)

    class Meta:
        db_table = 'PttTables'

class Ptttablesaccid(BaseModel):
    cog = TextField(db_column='COG', null=True)
    function = TextField(db_column='Function', null=True)
    position = IntegerField(db_column='Position', null=True)
    product = TextField(db_column='Product', null=True)
    synonymcode = TextField(db_column='SynonymCode', null=True)
    accid = CharField(db_column='accID', max_length=15, null=True)
    end = IntegerField(null=True)
    gene = TextField(null=True)
    gi = IntegerField(db_column='gi', primary_key=True)
    length = IntegerField(null=True)
    start = IntegerField(null=True)
    strand = CharField(max_length=1, null=True)

    class Meta:
        db_table = 'PttTablesAccID'

class Taxatable(BaseModel):
    orgclass = TextField(db_column='OrgClass', null=True)
    orgfamily = TextField(db_column='OrgFamily', null=True)
    orggenus = TextField(db_column='OrgGenus', null=True)
    orgname = TextField(db_column='OrgName', null=True)
    orgorder = TextField(db_column='OrgOrder', null=True)
    taxid = PrimaryKeyField()

    class Meta:
        db_table = 'TaxaTable'

class Taxidnpostable(BaseModel):
    position = IntegerField(null=True)
    taxid = IntegerField(null=True)

    class Meta:
        db_table = 'TaxidNPosTable'

class Taxidpostable(BaseModel):
    position = IntegerField(null=True)
    taxid = IntegerField(null=True)

    class Meta:
        db_table = 'TaxidPosTable'

class Temptable(BaseModel):
    cog = TextField(db_column='COG', null=True)
    function = TextField(db_column='Function', null=True)
    orgclass = TextField(db_column='OrgClass', null=True)
    orgfamily = TextField(db_column='OrgFamily', null=True)
    orggenus = TextField(db_column='OrgGenus', null=True)
    orgorder = TextField(db_column='OrgOrder', null=True)
    position = IntegerField(db_column='Position', null=True)
    product = TextField(db_column='Product', null=True)
    synonymcode = TextField(db_column='SynonymCode', null=True)
    accid = CharField(db_column='accID', max_length=15, null=True)
    end = IntegerField(null=True)
    gene = TextField(null=True)
    gi = IntegerField()
    homolog = TextField(null=True)
    length = IntegerField(null=True)
    method = TextField(null=True)
    organism = TextField(null=True)
    start = IntegerField(null=True)
    strand = CharField(max_length=1, null=True)

    class Meta:
        db_table = 'TempTable'

class Temporarytable(BaseModel):
    gi = IntegerField(null=True)

    class Meta:
        db_table = 'TemporaryTable'

class Gtahomologs(BaseModel):
    cog = TextField(db_column='COG', null=True)
    function = TextField(db_column='Function', null=True)
    position = IntegerField(db_column='Position', null=True)
    product = TextField(db_column='Product', null=True)
    synonymcode = TextField(db_column='SynonymCode', null=True)
    accid = CharField(db_column='accID', max_length=15, null=True)
    end = IntegerField(null=True)
    gene = TextField(null=True)
    gi = IntegerField()
    homolog = TextField(null=True)
    length = IntegerField(null=True)
    method = TextField(null=True)
    organism = TextField(null=True)
    start = IntegerField(null=True)
    strand = CharField(max_length=1, null=True)

    class Meta:
        db_table = 'gtahomologs'

class LargeCluster(BaseModel):
    accid = CharField(db_column='accID', max_length=15, null=True)
    cluster = IntegerField(null=True)
    cluster_count = TextField(null=True)
    organism = TextField(null=True)
    status = TextField(null=True)

    class Meta:
        db_table = 'large_cluster'

class TblastnHits(BaseModel):
    accid = CharField(db_column='accID', max_length=15, null=True)
    end = IntegerField(null=True)
    homolog = TextField(null=True)
    method = TextField(null=True)
    start = IntegerField(null=True)
    strand = TextField(null=True)

    class Meta:
        db_table = 'tblastn_hits'

