import fs from "fs";

export function readJson(path: string): any {
  return JSON.parse(fs.readFileSync(path, "utf8"));
}
