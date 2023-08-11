import fs from "fs";

export function readJson(path: string): any {
  const data = fs.readFileSync(path);
  return JSON.parse(data.toString("utf8"));
}
