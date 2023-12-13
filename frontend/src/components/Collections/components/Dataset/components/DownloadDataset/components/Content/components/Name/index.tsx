import { FormControl, FormLabel } from "@mui/material";
import { FC } from "react";

interface Props {
  name: string;
}

const Name: FC<Props> = ({ name }) => (
  <FormControl>
    <FormLabel>Name</FormLabel>
    <div data-testid="download-asset-name">{name}</div>
  </FormControl>
);

export default Name;
