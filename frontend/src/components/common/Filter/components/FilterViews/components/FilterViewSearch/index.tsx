import { IconNames } from "@blueprintjs/icons";
import {
  OnUpdateSearchValueFn,
  SetSearchValueFn,
} from "src/components/common/Filter/common/entities";
import { ViewSearch } from "./style";

interface Props {
  onUpdateSearchValue: OnUpdateSearchValueFn;
  setSearchValue: SetSearchValueFn;
}

export default function FilterViewSearch({
  onUpdateSearchValue,
  setSearchValue,
}: Props): JSX.Element {
  return (
    <ViewSearch
      leftIcon={IconNames.SEARCH}
      onChange={(changeEvent) =>
        onUpdateSearchValue(changeEvent, setSearchValue)
      }
      placeholder="Search"
    />
  );
}
