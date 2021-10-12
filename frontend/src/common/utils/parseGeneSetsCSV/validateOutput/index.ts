import {
  COMMENT_SYMBOL,
  Gene,
  GeneSet,
  Output,
  Result,
  RowInfo,
} from "./common";

export function validateOutput({
  errorMessages,
  databaseGeneSetNames,
  resolve,
  reject,
}: {
  databaseGeneSetNames: string[];
  errorMessages: Set<string>;
  resolve: null | ((value?: unknown) => void);
  reject: null | ((reason?: unknown) => void);
}) {
  return (parsed: Papa.ParseResult<string[]>) => {
    try {
      const { data, errors } = parsed;

      const output = buildOutput(data, errorMessages, databaseGeneSetNames);
      injectCSVParsingErrors(errors, errorMessages);

      resolve && resolve(output);
    } catch (error) {
      reject && reject(error);
    }
  };
}

function buildOutput(
  rows: string[][],
  errorMessages: Set<string>,
  databaseGeneSetNames: string[]
) {
  if (!rows.length) {
    return errorMessages.add("Seems to be an empty file");
  }

  const result: Result = new Map();

  let headerRow = null;
  let lastGeneSetName = "";

  const validators = [
    hasGeneSetDescription,
    hasDuplicateGeneSymbol,
    hasDuplicateGeneSetName,
    createHasDuplicateGeneSetNameInDatabase(databaseGeneSetNames),
  ];

  for (const [index, row] of Object.entries(rows)) {
    const rowInfo: RowInfo = {
      errorMessages,
      headerRow,
      lastGeneSetName,
      result,
      row,
      rowIndex: Number(index) + 1,
    };

    if (!headerRow) {
      headerRow = getHeaderRow(rowInfo);
      continue;
    }

    if (shouldSkip(rowInfo)) {
      continue;
    }

    for (const validator of validators) {
      const errorMessage = validator(rowInfo);

      if (errorMessage) {
        errorMessages.add(errorMessage);
      }
    }

    addResult(rowInfo);

    lastGeneSetName = row[0];
  }

  if (!headerRow) {
    errorMessages.add("Must provide the header row and in the required order");
  }

  const output: Output = Array.from(result.values()).map((geneSet) => {
    const genes = Array.from(geneSet.genes.values());

    return {
      ...geneSet,
      genes,
    };
  });

  return { errorMessages: Array.from(errorMessages), output };
}

function injectCSVParsingErrors(
  errors: Papa.ParseError[],
  errorMessages: Set<string>
) {
  if (!errors.length) return;

  for (const error of errors) {
    const { row, message } = error;

    errorMessages.add(`row: ${row}: ${message}`);
  }
}

const HEADER_ROW_KEYS = {
  GENE_DESCRIPTION: "gene_description",
  GENE_SET_DESCRIPTION: "gene_set_description",
  GENE_SET_NAME: "gene_set_name",
  GENE_SYMBOL: "gene_symbol",
};

function getHeaderRow({ row }: RowInfo) {
  const {GENE_DESCRIPTION, GENE_SET_DESCRIPTION, GENE_SET_NAME, GENE_SYMBOL} = HEADER_ROW_KEYS;

  const [gene_set_name, gene_set_description, gene_symbol, gene_description] = row;

  if (
    gene_set_name === GENE_SET_NAME &&
    gene_set_description === GENE_SET_DESCRIPTION &&
    gene_symbol === GENE_SYMBOL &&
    gene_description === GENE_DESCRIPTION
  ) {
    return row;
  } else {
    return null;
  }
}

const shouldSkipValidators = [
  hasGeneSetName,
  hasGeneSymbol,
  hasIllegalWhiteSpace,
];

function shouldSkip(rowInfo: RowInfo) {
  const { errorMessages } = rowInfo;

  if (isComment(rowInfo)) return true;

  let hasError = false;

  for (const validator of shouldSkipValidators) {
    const errorMessage = validator(rowInfo);

    if (errorMessage) {
      hasError = true;
      errorMessages.add(errorMessage);
    }
  }

  return hasError;
}

function isComment({ row }: RowInfo) {
  return row[0].trim()[0] === COMMENT_SYMBOL;
}

const ILLEGAL_WHITE_SPACE_REGEX = /^\s| {2}|[\v\t\r\n]|\s$/;

function hasIllegalWhiteSpace({ row, rowIndex }: RowInfo) {
  const isBad = row.some((value) => {
    return Boolean(value.match(ILLEGAL_WHITE_SPACE_REGEX));
  });

  return isBad
    ? `Row ${rowIndex}: The following row contains at least one illegal whitespace: Row "${row}".`
    : undefined;
}

function hasGeneSetName({ row, rowIndex }: RowInfo) {
  const [gene_set_name] = row;

  return gene_set_name ? undefined : `Row ${rowIndex}: Missing Gene Set Name`;
}

function hasGeneSetDescription({ row, rowIndex, result }: RowInfo) {
  const [gene_set_name, gene_set_description] = row;

  // (thuang): Ignore description after the first instance
  if (result.has(gene_set_name)) {
    return;
  }

  return gene_set_description
    ? undefined
    : `Row ${rowIndex}: Gene set: '${gene_set_name}' is missing Gene Set Description`;
}

function hasGeneSymbol({ row, rowIndex }: RowInfo) {
  const [gene_set_name, , gene_symbol] = row;

  return gene_symbol
    ? undefined
    : `Row ${rowIndex}: Gene set: '${gene_set_name}' is missing Gene Symbol`;
}

function hasDuplicateGeneSetName({
  lastGeneSetName,
  result,
  row,
  rowIndex,
}: RowInfo) {
  const [gene_set_name] = row;

  if (!lastGeneSetName || gene_set_name === lastGeneSetName) return;

  return result.has(gene_set_name)
    ? `Row ${rowIndex}: Gene set: '${gene_set_name}' was already used in this file. Please group all genes in a gene set together.`
    : undefined;
}

function createHasDuplicateGeneSetNameInDatabase(
  databaseGeneSetNames: string[]
) {
  return ({ row, rowIndex }: RowInfo) => {
    const [gene_set_name] = row;

    return databaseGeneSetNames.includes(gene_set_name)
      ? `Row ${rowIndex}: Gene set name: '${gene_set_name}' already exists in the database`
      : undefined;
  };
}

function hasDuplicateGeneSymbol({ result, row, rowIndex }: RowInfo) {
  const [gene_set_name, , gene_symbol] = row;

  const geneSet = result.get(gene_set_name);

  if (!geneSet) return;

  return geneSet.genes.has(gene_symbol)
    ? `Row ${rowIndex}: Gene symbol: '${gene_symbol}' already exists in the gene set: '${gene_set_name}'`
    : undefined;
}

function addResult({ row, result, headerRow }: RowInfo) {
  if (!headerRow) return;

  const [
    gene_set_name,
    gene_set_description,
    gene_symbol,
    gene_description,
    ...rest
  ] = row;

  const geneSet: GeneSet = result.get(gene_set_name) || {
    gene_set_description,
    gene_set_name,
    genes: new Map(),
  };

  const additional_params: Gene["additional_params"] = {};

  for (let i = 0; i < rest.length; i++) {
    const header = headerRow[i + 3];
    const value = rest[i];

    if (value) {
      additional_params[header] = value;
    }
  }

  const gene: Gene = {
    additional_params,
    gene_description,
    gene_symbol,
  };

  geneSet.genes.set(gene_symbol, gene);

  result.set(gene_set_name, geneSet);
}

export default validateOutput;
