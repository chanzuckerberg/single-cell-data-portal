import { DefaultMenuSelectOption } from "czifui";

export interface Organism {
  id: string;
  name: string;
}

export interface Filters {
  datasets?: DefaultMenuSelectOption[];
  developmentStages?: DefaultMenuSelectOption[];
  diseases?: DefaultMenuSelectOption[];
  ethnicities?: DefaultMenuSelectOption[];
  sexes?: DefaultMenuSelectOption[];
  tissues?: DefaultMenuSelectOption[];
}
