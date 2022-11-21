import { IconNames } from "@blueprintjs/icons";
import { SetSearchValueFn } from "src/components/common/Filter/components/FilterSearch/common/useFilterSearch";
import { ClearIconButton, ViewSearch } from "./style";
import { Icon } from "czifui";

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
      autoFocus
      className={className}
      leftIcon={IconNames.SEARCH}
      onChange={(changeEvent) => setSearchValue(changeEvent.target.value)}
      placeholder="Search"
      rightElement={
        searchValue ? (
          <ClearIconButton
            onClick={() => setSearchValue("")}
            sdsType="secondary"
          >
            <Icon
              sdsIcon="xMarkCircle"
              sdsSize="s" // maximum available IconNameToSize for xMarkCircle icon.
              sdsType="iconButton"
            />
          </ClearIconButton>
        ) : undefined
      }
      value={searchValue}
    />
  );
}
