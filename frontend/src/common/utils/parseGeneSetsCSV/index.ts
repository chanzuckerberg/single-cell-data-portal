import jschardet from "jschardet";
import Papa from "papaparse";
import { GeneSet as GeneSetEntity } from "src/common/entities";
import validateOutput from "./validateOutput";
import { COMMENT_SYMBOL } from "./validateOutput/common";

// (thuang): Please dynamic import this function, since it's heavy!
export async function parseGeneSetsCSV(
  file: File,
  databaseGeneSets: GeneSetEntity[] = []
) {
  const errorMessages = new Set<string>();
  const databaseGeneSetNames = getGeneSetNames(databaseGeneSets);

  // (thuang): Turn callback functions into promises
  let checkEncodingResolve: null | ((result?: unknown) => void) = null;
  let checkEncodingReject: null | ((reason?: unknown) => void) = null;
  let validateOutputResolve: null | ((result?: unknown) => void) = null;
  let validateOutputReject: null | ((reason?: unknown) => void) = null;

  const checkEncodingPromise = new Promise((resolve, reject) => {
    checkEncodingResolve = resolve;
    checkEncodingReject = reject;
  });

  const validateOutputPromise = new Promise((resolve, reject) => {
    validateOutputResolve = resolve;
    validateOutputReject = reject;
  });

  const fileReader = new FileReader();

  fileReader.onload = createCheckEncoding({
    errorMessages,
    reject: checkEncodingReject,
    resolve: checkEncodingResolve,
  });

  fileReader.readAsBinaryString(file);

  Papa.parse(file, {
    comments: COMMENT_SYMBOL,
    complete: validateOutput({
      databaseGeneSetNames,
      errorMessages,
      reject: validateOutputReject,
      resolve: validateOutputResolve,
    }),
  });

  const [result] = await Promise.all([
    validateOutputPromise,
    checkEncodingPromise,
  ]);

  return result;
}

function createCheckEncoding({
  errorMessages,
  reject,
  resolve,
}: {
  errorMessages: Set<string>;
  reject: null | ((reason?: unknown) => void);
  resolve: null | ((result?: unknown) => void);
}) {
  return (event: ProgressEvent<FileReader>) => {
    try {
      const result = String(event?.target?.result);

      const csvResult = result.split(/\r|\n|\r\n/);
      const encoding = jschardet.detect(csvResult.toString()).encoding;

      // (thuang): We only accept UTF and ASCII, since ASCII and UTF-8 are the
      // same for Latin letters
      // https://www.azavea.com/blog/2014/03/24/solving-unicode-problems-in-python-2-7/#1_str_is_for_bytes_not_strings
      if (!encoding?.includes("UTF") && !encoding?.includes("ascii")) {
        if (encoding) {
          errorMessages.add(
            `The CSV encoding is '${encoding}', instead of a 'UTF'`
          );
        } else {
          errorMessages.add(
            "Could not detect the encoding. Please choose a UTF"
          );
        }
      }

      resolve && resolve();
    } catch (error) {
      reject && reject(error);
    }
  };
}

function getGeneSetNames(geneSets: GeneSetEntity[]): string[] {
  return geneSets.map((geneSet) => geneSet.name);
}

export default parseGeneSetsCSV;
