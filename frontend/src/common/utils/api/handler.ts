import type { NextApiRequest, NextApiResponse } from "next";
import readJson from "./readJson";

const getHandler = (path: string) => {
  const data = readJson(path);

  async function handler(req: NextApiRequest, res: NextApiResponse) {
    try {
      const cellTypeId = req.query.cellTypeId;
      if (data[cellTypeId as keyof typeof data])
        res.status(200).json(data[cellTypeId as keyof typeof data]);
      else res.status(404);
    } catch (error) {
      console.error(`Error fetching data: ${error}`);
      res.status(500).json({ error: "Error fetching data" });
    }
  }
  return handler;
};

export default getHandler;
