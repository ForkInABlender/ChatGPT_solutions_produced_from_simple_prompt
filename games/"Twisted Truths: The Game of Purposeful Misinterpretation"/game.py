# Dylan Kenneth Eliot & GPT-4o ( Alpha Edition )

"""

### Twisted Truths: The Game of Purposeful Misinterpretation

#### Objective
The objective of the game is to accumulate the most points by creatively misinterpreting statements and convincing other players of your twisted truths.

#### Players
3-6 players

#### Setup
1. Each player is given a set of 5 Misinterpretation Cards.
2. A deck of Statement Cards is shuffled and placed in the center of the table.
3. Each player draws 3 Statement Cards to start the game.

#### Gameplay
1. **Starting the Round:** The player who most recently told a funny story starts the game as the Reader.
2. **Reading a Statement:** The Reader selects a Statement Card from their hand and reads it aloud to the group.
3. **Misinterpretation Phase:** Each player, including the Reader, selects a Misinterpretation Card from their hand and writes a creative and plausible
                                 misinterpretation of the Statement on it.
4. **Presenting Misinterpretations:** Players take turns presenting their misinterpretations to the group.
5. **Voting:** After all misinterpretations are presented, players vote on which misinterpretation they find the most convincing or entertaining. Players
                cannot vote for their own misinterpretation.
6. **Scoring:**
   - The player whose misinterpretation receives the most votes earns 3 points.
   - If there is a tie, each tied player earns 2 points.
   - The Reader earns 1 point if no one guesses the correct interpretation.

#### End of the Round
1. After scoring, each player draws one new Statement Card and one new Misinterpretation Card.
2. The next player clockwise becomes the Reader, and the new round begins.

#### Winning the Game
The game continues for a predetermined number of rounds (e.g., 10 rounds) or until the Statement Card deck is exhausted. The player with the most points at
 the end of the game wins.

#### Misinterpretation Cards Examples
1. "Change the context to a different time period."
2. "Swap the subject and the object."
3. "Add a humorous twist."
4. "Incorporate a famous quote or saying."
5. "Change the setting to a fantasy world."

Players are encouraged to be as creative and convincing as possible while misinterpreting the statements to make the game more engaging and fun.



"""


import random

class Player:
    def __init__(self, name):
        self.name = name
        self.statement_cards = []
        self.misinterpretation_cards = []
        self.points = 0

    def draw_statement_card(self, deck):
        self.statement_cards.append(deck.draw_card())

    def draw_misinterpretation_card(self, deck):
        self.misinterpretation_cards.append(deck.draw_card())

    def play_statement_card(self):
        return self.statement_cards.pop(0)

    def play_misinterpretation_card(self):
        return self.misinterpretation_cards.pop(0)

    def __str__(self):
        return self.name

class Deck:
    def __init__(self, cards):
        self.cards = cards
        random.shuffle(self.cards)

    def draw_card(self):
        return self.cards.pop(0) if self.cards else None

def setup_game():
    statement_cards = [
        "The sky is blue.", "Dogs are better than cats.", 
        "Chocolate is better than vanilla.", "Summer is the best season.",
        "Books are better than movies."
    ]
    misinterpretation_cards = [
        "Change the context to a different time period.", "Swap the subject and the object.",
        "Add a humorous twist.", "Incorporate a famous quote or saying.",
        "Change the setting to a fantasy world."
    ]

    statement_deck = Deck(statement_cards)
    misinterpretation_deck = Deck(misinterpretation_cards * 6)  # Assume plenty of misinterpretation cards

    players = [Player(f"Player {i + 1}") for i in range(3)]
    for player in players:
        for _ in range(3):
            player.draw_statement_card(statement_deck)
        for _ in range(5):
            player.draw_misinterpretation_card(misinterpretation_deck)
    
    return players, statement_deck, misinterpretation_deck

def play_round(players, statement_deck, misinterpretation_deck):
    reader = players.pop(0)
    statement = reader.play_statement_card()
    print(f"{reader} reads: {statement}")
    
    misinterpretations = {}
    for player in players:
        misinterpretation = player.play_misinterpretation_card()
        print(f"{player} misinterprets: {misinterpretation}")
        misinterpretations[player] = misinterpretation

    # Voting phase (simplified to random choice for demonstration)
    votes = {player: 0 for player in players}
    for player in players:
        vote = random.choice(list(votes.keys()))
        if vote != player:
            votes[vote] += 1

    # Scoring phase
    max_votes = max(votes.values())
    for player, vote_count in votes.items():
        if vote_count == max_votes:
            player.points += 3
    
    reader.points += 1 if max_votes == 0 else 0

    # Draw new cards
    reader.draw_statement_card(statement_deck)
    for player in players:
        player.draw_misinterpretation_card(misinterpretation_deck)
    
    players.append(reader)  # Reader goes to the end of the list

def main():
    players, statement_deck, misinterpretation_deck = setup_game()
    
    rounds = 10
    for _ in range(rounds):
        play_round(players, statement_deck, misinterpretation_deck)

    print("\nFinal Scores:")
    for player in players:
        print(f"{player}: {player.points} points")

if __name__ == "__main__":
    main()
