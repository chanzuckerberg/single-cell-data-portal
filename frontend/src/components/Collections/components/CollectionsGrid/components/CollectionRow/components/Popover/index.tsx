import {
  Button,
  Intent,
  Popover as PopoverRaw,
  Position,
} from "@blueprintjs/core";
import React, { FC } from "react";
import { LeftAlignedDetailsCell } from "../common/style";
import { ContentWrapper, FieldValues } from "./style";

interface Props {
  values: string[];
}

const Popover: FC<Props> = ({ values }) => {
  return (
    <LeftAlignedDetailsCell>
      <FieldValues>
        {values[0]}
        <br />
        {values[1]}
      </FieldValues>
      {values.length > 2 && (
        <PopoverRaw
          position={Position.BOTTOM}
          content={
            <ContentWrapper>
              {values.map((val, idx) => (
                <FieldValues key={val}>
                  {val}
                  {idx !== values.length - 1 && <br />}
                </FieldValues>
              ))}
            </ContentWrapper>
          }
        >
          <Button small outlined intent={Intent.PRIMARY}>
            +{values.length - 2}
          </Button>
        </PopoverRaw>
      )}
    </LeftAlignedDetailsCell>
  );
};

export default Popover;
