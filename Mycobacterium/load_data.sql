USE mycobacterium;

DROP TABLE samples, dst, corres, sras, resistance;

CREATE TABLE samples(Specimen varchar(255) NOT NULL,Patient varchar(255) NOT NULL, Bio_Project varchar(255)  NOT NULL,Material  varchar(255) NOT NULL,Specimen_Collected_Date Date NOT NULL,Bio_Sample varchar(255) NOT NULL,Lineage varchar(255) NOT NULL,Octal_Spoligotype varchar(255) NOT NULL, PRIMARY KEY (Specimen));

CREATE TABLE sras(Specimen varchar(255) NOT NULL, Sequence_Read_Archive varchar(255) NOT NULL, PRIMARY KEY(Specimen, Sequence_Read_Archive));

CREATE TABLE dst(Specimen varchar(255) NOT NULL, test varchar(255), Date date, H varchar(255) NOT NULL, R varchar(255) NOT NULL, S varchar(255) NOT NULL, E varchar(255) NOT NULL, Ofx varchar(255) NOT NULL, Cm varchar(255) NOT NULL, Am varchar(255) NOT NULL, Km varchar(255) NOT NULL, Z varchar(255) NOT NULL, Lfx varchar(255) NOT NULL, Mfx varchar(255) NOT NULL, Pas varchar(255) NOT NULL, Pto varchar(255) NOT NULL, Cs varchar(255) NOT NULL, Amx_Clv varchar(255) NOT NULL, Mb varchar(255) NOT NULL, Dld varchar(255) NOT NULL, Bdq varchar(255) NOT NULL, Ipm_Cln varchar(255) NOT NULL, Lzd varchar(255) NOT NULL, Cfz varchar(255) NOT NULL, Clr varchar(255) NOT NULL, Ft varchar(255) NOT NULL, AG_CP varchar(255) NOT NULL, Action varchar(255) NOT NULL, PRIMARY KEY (Specimen, test));

CREATE TABLE corres(Specimen varchar(255) NOT NULL, file varchar(255), PRIMARY KEY(Specimen));

CREATE TABLE resistance(antibiotic varchar(255) NOT NULL, gene varchar(255), PRIMARY KEY(antibiotic, gene));

LOAD DATA LOCAL INFILE 'samples.csv' INTO TABLE samples FIELDS TERMINATED BY '\t' ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS;
LOAD DATA LOCAL INFILE 'test.csv' INTO TABLE dst FIELDS TERMINATED BY '\t' ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS;
LOAD DATA LOCAL INFILE 'file_correspondance.csv' INTO TABLE corres FIELDS TERMINATED BY '\t' ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS;
LOAD DATA LOCAL INFILE 'sras.csv' INTO TABLE sras FIELDS TERMINATED BY '\t' ENCLOSED BY '"' LINES TERMINATED BY '\n' IGNORE 1 ROWS;
LOAD DATA LOCAL INFILE '~/Documents/scripts/Mycobacterium/resistance_genes.csv' INTO table resistance FIELDS TERMINATED BY '\t' ENCLOSED BY '"' LINES TERMINATED BY '\n';

