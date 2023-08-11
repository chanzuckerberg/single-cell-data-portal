import fs from "fs";

export function readJson(path: string): unknown {
  const data = fs.readFileSync(path);
  return JSON.parse(data.toString("utf8"));
}
