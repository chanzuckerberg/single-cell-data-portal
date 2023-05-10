import type { NextApiRequest, NextApiResponse } from "next";
import fs from "fs";

let allSourceData = {};
fs.readFile(
  "src/views/CellCards/common/sourceDataFixture.json",
  "utf8",
  (err, jsonString) => {
    if (err) {
      return;
    }
    try {
      allSourceData = JSON.parse(jsonString);
    } catch (err) {
      console.log("Error parsing JSON string:", err);
    }
  }
);

async function handler(req: NextApiRequest, res: NextApiResponse) {
  try {
    const cellTypeId = req.query.cellTypeId as string;
    if (allSourceData[cellTypeId as keyof typeof allSourceData])
      res
        .status(200)
        .json(allSourceData[cellTypeId as keyof typeof allSourceData]);
    else res.status(500).json({ error: "Source data not available yet." });
  } catch (error) {
    console.error(`Error fetching source data info: ${error}`);
    res.status(500).json({ error: "Error fetching source data info" });
  }
}

export default handler;
