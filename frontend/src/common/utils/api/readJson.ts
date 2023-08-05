import fs from "fs";
import zlib from "zlib";

export function readJson(path: string): any {
  const data = fs.readFileSync(path);

  // Check if the file is gzip compressed by examining the first two bytes (magic numbers)
  const isGzipCompressed =
    (data[0] === 0x1f && data[1] === 0x8b) || path.endsWith(".gz");

  let jsonData: string;

  if (isGzipCompressed) {
    const decompressedData = zlib.gunzipSync(data);
    jsonData = decompressedData.toString("utf8");
  } else {
    jsonData = data.toString("utf8");
  }

  return JSON.parse(jsonData);
}
