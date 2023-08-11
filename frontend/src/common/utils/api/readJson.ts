import fs from "fs";
import AdmZip from "adm-zip";

// handle unzipping of cellguide fixtures
const CELLGUIDE_FIXUTRES_ZIPPED_FILE = `src/views/CellGuide/common/fixtures/cellguide_data.zip`;

function unzipFolder(zipFilePath: string, outputFolderPath: string): void {
  const zip = new AdmZip(zipFilePath);
  zip.extractAllTo(outputFolderPath, true); // true for overwrite
}
unzipFolder(CELLGUIDE_FIXUTRES_ZIPPED_FILE, ".");

export function readJson(path: string): any {
  const data = fs.readFileSync(path);
  return JSON.parse(data.toString("utf8"));
}
