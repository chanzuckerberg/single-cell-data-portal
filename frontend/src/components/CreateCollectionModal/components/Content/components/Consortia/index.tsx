import { useState } from "react";
import Dropdown, {
  Value as ConsortiaValue,
} from "src/components/common/Form/Dropdown";

// TODO(cc) grab BE consortia values.
const CONSORTIA_OPTIONS = [
  {
    name: "Allen Institute for Brain Science",
  },
  {
    name: "BRAIN Initiative",
  },
  {
    name: "CZ Biohub",
  },
  {
    name: "CZI Single-Cell Biology",
  },
  {
    name: "European Union’s Horizon 2020",
  },
  {
    name: "GenitoUrinary Development Molecular Anatomy Project (GUDMAP)",
  },
  {
    name: "Human Tumor Atlas Network (HTAN)",
  },
];

export default function Consortia(): JSX.Element {
  const [consortia, setConsortia] = useState<ConsortiaValue>(null);

  const consortiaOnChange = (newConsortia: ConsortiaValue) => {
    setConsortia(newConsortia);
  };

  return (
    <Dropdown
      label="Select Consortia"
      multiple
      onChange={consortiaOnChange}
      optionalField
      options={CONSORTIA_OPTIONS}
      text="Consortia"
      value={consortia}
    />
  );
}
