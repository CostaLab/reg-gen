.header on
.mode csv

--TFs with no entry in MATRIX_ANNOTATION or TAX will be excluded

SELECT
  trim(replace((m.BASE_ID || '.' || m.VERSION || '.' || m.NAME),'/','_')) AS TfName,
  trim(replace(m.NAME, '::', '+')) AS GeneName,
  trim(replace(GROUP_CONCAT(DISTINCT CASE WHEN a.TAG = 'family' THEN a.VAL END),',',';')) AS Family,
  trim(replace(GROUP_CONCAT(DISTINCT p.ACC),',',';')) AS UniProt,
  trim(replace(GROUP_CONCAT(DISTINCT CASE WHEN a.tag = 'type' THEN a.val END),',',';')) AS Source,
  trim(GROUP_CONCAT(DISTINCT CASE WHEN a.tag = 'tax_group' THEN a.val END)) AS TaxGroup,
  trim(s.SPECIES) AS Species
FROM MATRIX m
  JOIN (SELECT *
        FROM MATRIX_ANNOTATION a2
        ORDER BY ID) a
    ON m.ID = a.ID
  LEFT OUTER JOIN MATRIX_PROTEIN p
    ON m.ID = p.ID
  LEFT OUTER JOIN (SELECT *
        FROM MATRIX_SPECIES sp
          JOIN TAX t
            ON sp.TAX_ID = t.TAX_ID
        ORDER BY sp.ID) s
    ON m.ID = s.ID
WHERE m.COLLECTION = 'CORE'
GROUP BY m.ID, m.COLLECTION, m.BASE_ID, m.VERSION, m.NAME
ORDER BY m.ID;
