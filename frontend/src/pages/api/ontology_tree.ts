import type { NextApiRequest, NextApiResponse } from "next";
import ontologyTree from "src/views/CellGuide/common/fixtures/ontologyTree.json";

async function handler(_: NextApiRequest, res: NextApiResponse) {
  try {
    res.status(200).json(ontologyTree);
  } catch (error) {
    console.error(`Error fetching data: ${error}`);
    res.status(500).json({ error: "Error fetching data" });
  }
}

export default handler;
