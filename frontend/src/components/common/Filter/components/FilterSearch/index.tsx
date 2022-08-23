import { IconNames } from "@blueprintjs/icons";
import { SetSearchValueFn } from "src/components/common/Filter/common/entities";
import { ViewSearch } from "./style";

interface Props {
  className?: string;
  searchValue: string;
  setSearchValue: SetSearchValueFn;
}

export default function FilterSearch({
  className,
  searchValue,
  setSearchValue,
}: Props): JSX.Element {
  return (
    <ViewSearch
      className={className}
      leftIcon={IconNames.SEARCH}
      onChange={(changeEvent) => setSearchValue(changeEvent.target.value)}
      placeholder="Search"
      value={searchValue}
    />
  );
}
